# src/Analysis.jl
 
########################################################################
#
#  Analysis.jl
#
#  Analyse-functies voor de FisherGeometrics informatiefunctionaal:
#
#    I(R) = S_F(ρ(R)) / D²_Bures(ρ(R), ρ*)
#
#  Het minimum van I(R) correspondeert met de compactificatieschaal
#  die het hiërarchiprobleem oplost (conjectuur: R* ≈ 4260).
#
########################################################################
 
using LinearAlgebra, Printf
 
# ── Informatiefunctionaal via exacte Christoffel ──────────────────
 
"""
    compute_I_R(R_values; n=6, verbose=true) -> Vector{NamedTuple}
 
Berekent de Fisher-informatie functionaal I(R) = S_F(ρ(R)) / D²(ρ(R),ρ*)
voor een reeks R-waarden via exacte metriekafgeleiden (`christoffel_rotate`).
 
De S_F berekening gebruikt alleen de ΓΓ-bijdrage aan de Riemann-tensor.
Dit is exact bij ρ* (Theorem E: ∂Γ = 0) en een goede benadering voor R > 200.
Voor kleine R (< 75) worden de resultaten onbetrouwbaar.
 
# Argumenten
- `R_values`: vector van compactificatiestralen
- `n`: dimensie van de toestandsruimte (standaard 6 voor su(6))
 
# Geeft terug
Vector van NamedTuples met velden `(R, D2, SF, I)`.
 
# Gebruik
```julia
T       = su_basis(6)
results = compute_I_R([100, 200, 500, 1000, 4260])
print_I_R_summary(results)
```
 
# Zie ook
- `compute_I_R_exact`: oudere versie via SF_exact (langzamer)
- `compute_I_R_path`: pad-evolutie met entropie
"""
function compute_I_R(R_values::Vector{<:Real}; n::Int=6,
                     verbose::Bool=true)
    T = su_basis(n)
    results = NamedTuple[]
 
    verbose && println("═"^60)
    verbose && println("I(R) = S_F(ρ(R)) / D²   [via christoffel_rotate]")
    verbose && println("═"^60)
    verbose && @printf("  %8s  %10s  %8s  %12s\n", "R", "S_F", "D²", "I(R)")
    verbose && println("  " * "─"^42)
 
    for R in R_values
        p  = rho_KK_eigenvalues(R; n=n)
        D2 = D2_bures_KK(R; n=n)
 
        # S_F via exacte ΓΓ-Riemann
        Γ  = christoffel_rotate(p, T)
        N  = length(T)
 
        # Metriek inverse
        G  = zeros(Float64, N, N)
        for a in 1:N, b in 1:N
            G[a,b] = sum(real(conj(T[a][i,j])*T[b][i,j])/(2*(p[i]+p[j]))
                         for i in 1:n, j in 1:n)
        end
        Gi = pinv(G; atol=1e-10)
 
        # Riemann ΓΓ → Ricci → S_F
        R_GG = zeros(Float64, N, N, N, N)
        for e in 1:N, a in 1:N, b in 1:N, c in 1:N
            R_GG[e,a,b,c] = sum(Γ[f,b,c]*Γ[e,a,f] - Γ[f,a,c]*Γ[e,b,f]
                                 for f in 1:N)
        end
        Ric_GG = [sum(R_GG[b,a,b,c] for b in 1:N) for a in 1:N, c in 1:N]
        SF = sum(Gi[a,c]*Ric_GG[a,c] for a in 1:N, c in 1:N)
 
        I_R = D2 > 1e-12 ? SF / D2 : Inf
        push!(results, (R=R, D2=D2, SF=SF, I=I_R))
 
        verbose && @printf("  %8.0f  %+10.2f  %8.5f  %12.1f\n",
                           R, SF, D2, I_R)
    end
 
    if verbose && !isempty(results)
        idx = argmin([r.I for r in results])
        println("  " * "─"^42)
        @printf("  Minimum bij R* ≈ %.0f  (I = %.1f)\n",
                results[idx].R, results[idx].I)
    end
    return results
end
 
"""
    compute_I_R_exact(R_values; n=6) -> Vector{NamedTuple}
 
Oudere versie van `compute_I_R` via `SF_exact`.
Bewaard voor achterwaartse compatibiliteit.
"""
function compute_I_R_exact(R_values::Vector{<:Real}; n::Int=6)
    T = su_basis(n)
    results = NamedTuple[]
    for R in R_values
        ρ  = rho_KK(R; n=n)
        D2 = D2_bures_KK(R; n=n)
        SF = SF_exact(ρ, T)
        I_R = D2 > 1e-12 ? SF / D2 : Inf
        push!(results, (R=R, D2=D2, SF=SF, I=I_R))
    end
    return results
end
 
# ── Output ────────────────────────────────────────────────────────
 
"""
    print_I_R_summary(results)
 
Drukt een geformatteerde tabel af van `compute_I_R` resultaten.
"""
function print_I_R_summary(results::Vector{<:NamedTuple})
    println("═"^60)
    println("Fisher Informatie Analyse: I(R) = S_F(ρ(R)) / D²")
    println("═"^60)
    for r in results
        @printf("R = %8.1f | S_F = %+9.3f | D² = %8.5f | I = %10.1f\n",
                r.R, r.SF, r.D2, r.I)
    end
    I_vals = [r.I for r in results]
    idx = argmin(I_vals)
    println("─"^60)
    @printf("Minimum I(R) bij R* ≈ %.1f  (I = %.1f)\n",
            results[idx].R, results[idx].I)
end
 
# ── Pad-analyse ───────────────────────────────────────────────────
 
"""
    compute_I_R_path(ρ_start, direction, steps, step_size, T)
        -> Vector{Tuple{Float64,Float64,Float64}}
 
Berekent S_F (ΓΓ) en Von Neumann entropie langs een lineair pad
in de toestandsruimte D₆:
 
    ρ(t) = (ρ_start + t × direction) / Tr(...)
 
Geeft een vector van tuples `(t, S_F, S_VN)`.
 
# Gebruik
```julia
T      = su_basis(6)
ρ_mass = rho_KK(100.0)
ρ★     = Matrix{ComplexF64}(I,6,6)/6
results = compute_I_R_path(ρ_mass, ρ★ - ρ_mass, 10, 0.1, T)
for (t, sf, svn) in results
    @printf("t=%.1f  S_F=%.3f  S_VN=%.4f\\n", t, sf, svn)
end
```
"""
function compute_I_R_path(ρ_start::AbstractMatrix,
                           direction::AbstractMatrix,
                           steps::Int,
                           step_size::Float64,
                           T::Vector)
    N = length(T)
    n = size(ρ_start, 1)
    results = Tuple{Float64,Float64,Float64}[]

    for i in 0:steps
        t   = i * step_size
        ρ_c = Matrix{ComplexF64}(ρ_start + t * direction)
        ρ_c /= tr(ρ_c)

        # Eigenwaarden voor christoffel_rotate
        p = sort(max.(real(eigvals(Hermitian(ρ_c))), 1e-15))

        # S_F via exacte ΓΓ (christoffel_rotate — geen SF_GG)
        Γ  = christoffel_rotate(p, T)
        G  = zeros(Float64, N, N)
        for a in 1:N, b in 1:N
            G[a,b] = sum(real(conj(T[a][ii,jj])*T[b][ii,jj])/(2*(p[ii]+p[jj]))
                         for ii in 1:n, jj in 1:n)
        end
        Gi = pinv(G; atol=1e-10)

        R_GG = zeros(Float64, N, N, N, N)
        for e in 1:N, a in 1:N, b in 1:N, c in 1:N
            R_GG[e,a,b,c] = sum(Γ[f,b,c]*Γ[e,a,f] - Γ[f,a,c]*Γ[e,b,f]
                                 for f in 1:N)
        end
        Ric = [sum(R_GG[b,a,b,c] for b in 1:N) for a in 1:N, c in 1:N]
        sf  = sum(Gi[a,c]*Ric[a,c] for a in 1:N, c in 1:N)

        # Von Neumann entropie
        svn = -sum(v * log(v) for v in p if v > 1e-15)

        push!(results, (t, sf, svn))
        @printf("  t=%.2f  S_F=%+10.3f  S_VN=%.4f\n", t, sf, svn)
    end
    return results
end
 
# ── Berry-fase analyse ────────────────────────────────────────────
 
"""
    run_fase_test(R_range, β, M1, M2, J) -> Float64
 
Berekent de geometrische (Berry) fase langs een radiaal pad R_range
voor een KK-model met inverse temperatuur β en massa-parameters M1, M2, J.
 
Een niet-nul fase duidt op topologische lading of kromming.
"""
function run_fase_test(R_range, β, M1, M2, J)
    total_phase = 0.0
    for i in 1:length(R_range)-1
        R1, R2 = R_range[i], R_range[i+1]
 
        ρ1 = gibbs_state_expanded(M1, M2, J/R1, β)
        ρ2 = gibbs_state_expanded(M1, M2, J/R2, β)
 
        val1 = 0.5 * exp(im * (2π * R1 / 4260))
        val2 = 0.5 * exp(im * (2π * R2 / 4260))
        ρ1[1,2] += val1; ρ1[2,1] += conj(val1)
        ρ2[1,2] += val2; ρ2[2,1] += conj(val2)
        ρ1 /= tr(ρ1); ρ2 /= tr(ρ2)
 
        _, V1 = eigen(ρ1)
        _, V2 = eigen(ρ2)
 
        overlap = dot(V1[:,end], V2[:,end])
        total_phase += angle(overlap)
    end
    return total_phase
end
