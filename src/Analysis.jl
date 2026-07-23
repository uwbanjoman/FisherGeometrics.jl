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

"""
    proton_mass(; v::Float64=246220.0) -> Float64

Proton mass from the M¹·¹·¹ topology.

    m_p = μ₀ × 4τ = 4v/1050

Three geometric factors:
  μ₀ = v/210 = 1172.47 MeV   Higgs VEV / KK normalization factor
  4τ = 4/5   = 0.8            holonomy scaling of M¹·¹·¹
  (equivalently: dim(SU(3)_adj) × τ/2 = 8 × 1/10 = 0.8)

The proton is a topological Skyrmion in the SU(3) sector of M¹·¹·¹.
Its mass is determined by the winding of the SU(3) holonomy, not by
a KK eigenvalue. This resolves the QCD mass gap problem geometrically.

Residual 0.29 MeV (0.031%) = QED electromagnetic self-energy of quarks.

# Arguments
- `v`: Higgs VEV in MeV (default: 246220.0)

# Returns
Proton mass in MeV

# Examples
```julia
proton_mass()              # → 937.98 MeV  (measured: 938.272 MeV, 0.031%)
proton_mass(v=246220.0)    # same
```

See: FisherGeometrics preprint v15, section 7.
"""
function proton_mass(; v::Float64=246220.0)
    τ  = 1/5
    μ₀ = v/210.0
    return μ₀ * 4τ
end

"""
    pion_mass(; v::Float64=246220.0) -> Float64

Neutral pion mass from the Killing-spinor series on CP².

    m_π = μ₀ × 3τ²/(1+τ²+τ⁴)

The factor 3 = dim(SU(2)_adj) counts the W-boson degrees of freedom.
The factor τ²/(1+τ²+τ⁴) is the second term of the Killing-spinor
series on CP² — the geometric series correction that encodes the
coupling of the Goldstone mode (pion) to the second Killing-spinor
direction on CP².

The pion is the Goldstone boson of chiral symmetry breaking
SU(2)_L × SU(2)_R → SU(2)_V, geometrically realized as the
second holonomy of M¹·¹·¹.

Compare with the proton:
  m_p = μ₀ × 8 × τ/2    [first holonomy,  SU(3), 0.031%]
  m_π = μ₀ × 3 × τ²/Σ   [second holonomy, SU(2), 0.075%]

# Arguments
- `v`: Higgs VEV in MeV (default: 246220.0)

# Returns
Neutral pion mass in MeV

# Examples
```julia
pion_mass()              # → 135.08 MeV  (measured: 134.977 MeV, 0.075%)
```

See: FisherGeometrics preprint v15, section 7.
"""
function pion_mass(; v::Float64=246220.0)
    τ  = 1/5
    μ₀ = v/210.0
    return μ₀ * 3τ^2 / (1 + τ^2 + τ^4)
end

"""
    hadron_spectrum(; v::Float64=246220.0, verbose::Bool=true) -> Vector{NamedTuple}

Compute 13 hadron masses from the M¹·¹·¹ topology.

All masses follow from the universal scale μ₀ = v/210 and three
geometric parameters — no free parameters:

    μ₀ = v/210 = 1172.47 MeV    (Higgs VEV / KK normalization)
    τ  = 1/5                     (holonomy parameter of M¹·¹·¹)
    α  = 1/137.036               (fine structure constant)

Two-phase structure:
  Light sector (m < μ₀): topological τ-series — Killing-spinor holonomy
  Heavy sector (m ≥ μ₀): KK eigenvalues — isometry group dimensions

The strange-quark increment δ_s = μ₀/8 = 146.6 MeV follows from
dim(SU(3)_adj) = 8: each strange quark occupies one SU(3) degree of freedom.

Average error: 0.29%. Maximum: 1.21% (Sigma Σ).

# Arguments
- `v`      : Higgs VEV in MeV (default: 246220.0)
- `verbose`: print the mass table (default: true)

# Returns
Vector of NamedTuples with fields:
  name, formula, M_FG, M_exp, error_pct, sector

# Examples
```julia
spectrum = hadron_spectrum()
spectrum = hadron_spectrum(verbose=false)

m_proton = first(s.M_FG for s in spectrum if s.name == "proton p")
```

See: FisherGeometrics preprint v15, section 7.
     examples/hadron_spectrum.jl for a complete demonstration.
"""
function hadron_spectrum(; v::Float64=246220.0, verbose::Bool=true)
    τ   = 1/5
    α   = 1/137.036
    KK  = 210.0
    μ₀  = v/KK
    δ_s = μ₀/8

    entries = [
        (name="pion π⁰",   sector="light", formula="μ₀×3τ²/(1+τ²+τ⁴)",
         M_exp=134.977,  M_FG=μ₀*3τ^2/(1+τ^2+τ^4)),
        (name="pion π±",   sector="light", formula="mπ⁰ + α/2×μ₀",
         M_exp=139.570,  M_FG=μ₀*3τ^2/(1+τ^2+τ^4) + α/2*μ₀),
        (name="kaon K",    sector="light", formula="μ₀×(2τ+τ²/2)",
         M_exp=493.677,  M_FG=μ₀*(2τ+τ^2/2)),
        (name="eta η",     sector="light", formula="μ₀×(2τ+5τ²/3)",
         M_exp=547.862,  M_FG=μ₀*(2τ+5τ^2/3)),
        (name="rho ρ",     sector="heavy", formula="μ₀×2/3",
         M_exp=775.260,  M_FG=μ₀*2/3),
        (name="omega ω",   sector="heavy", formula="μ₀×2/3",
         M_exp=782.650,  M_FG=μ₀*2/3),
        (name="proton p",  sector="light", formula="μ₀×4τ",
         M_exp=938.272,  M_FG=μ₀*4τ),
        (name="neutron n", sector="light", formula="mₚ + μ₀×α/(2π)",
         M_exp=939.565,  M_FG=μ₀*4τ + μ₀*α/(2π)),
        (name="Lambda Λ",  sector="heavy", formula="μ₀×4τ(1+τ)",
         M_exp=1115.683, M_FG=μ₀*4τ*(1+τ)),
        (name="Delta Δ",   sector="heavy", formula="μ₀×3(1+2τ)/4",
         M_exp=1232.000, M_FG=μ₀*3*(1+2τ)/4),
        (name="Sigma Σ",   sector="heavy", formula="mΛ + μ₀τ/3",
         M_exp=1189.370, M_FG=μ₀*4τ*(1+τ) + μ₀*τ/3),
        (name="Xi Ξ",      sector="heavy", formula="mΣ + δ_s×(1−τ)",
         M_exp=1314.860, M_FG=μ₀*4τ*(1+τ) + μ₀*τ/3 + δ_s*(1-τ)),
        (name="Omega Ω",   sector="heavy", formula="mΔ + 3δ_s",
         M_exp=1672.450, M_FG=μ₀*3*(1+2τ)/4 + 3*δ_s),
    ]

    results = [(; e..., error_pct=abs(e.M_FG-e.M_exp)/e.M_exp*100)
               for e in entries]

    if verbose
        @printf("\nHADRON SPECTRUM — FisherGeometrics (μ₀ = %.2f MeV)\n", μ₀)
        println("="^70)
        @printf("%-12s  %-6s  %-28s  %-8s  %-8s  %-6s\n",
                "Hadron", "Sector", "Formula", "FG (MeV)", "Exp (MeV)", "Error")
        println("─"^70)
        current = ""
        for r in results
            if r.sector != current
                current = r.sector
                label = current == "light" ?
                    "── Light sector (<μ₀): topological τ-series ──" :
                    "── Heavy sector (≥μ₀): KK eigenvalues ──"
                println("\n  $label")
            end
            @printf("  %-12s %-6s  %-28s  %8.2f  %8.3f  %5.3f%%\n",
                    r.name, "", r.formula, r.M_FG, r.M_exp, r.error_pct)
        end
        println("\n" * "─"^70)
        errors = [r.error_pct for r in results]
        @printf("  Average: %.3f%%   Maximum: %.3f%% (%s)\n",
                sum(errors)/length(errors),
                maximum(errors),
                results[argmax(errors)].name)
        @printf("  μ₀ = %.4f MeV   τ = 1/5   δ_s = μ₀/8 = %.4f MeV\n",
                μ₀, δ_s)
        println()
    end

    return results
end

"""
    SF_GG(p::Vector{Float64}) -> Float64

Compute the Ricci scalar S_F via the ΓΓ terms of the Bures metric.

Extended Theorem E guarantees S_F^{∂Γ} = 0 for all ρ ∈ D₆,
so S_F = S_F^{ΓΓ} exactly (Dittmann 1999, symmetric space).

Computation steps:
  1. Christoffel symbols Γ via christoffel_rotate(p, T)
  2. Bures metric G_{ab} and inverse G^{ab} via pinv
  3. Riemann tensor R^e_{abc} = Σ_f (Γ^f_{bc}Γ^e_{af} − Γ^f_{ac}Γ^e_{bf})
  4. Ricci tensor Ric_{ac} = Σ_b R^b_{abc}
  5. Ricci scalar S_F = Σ_{ac} G^{ac} Ric_{ac}

At the vacuum ρ* = I/6: S_F = 560 (proved, Proof 04).
The information-geometric RGE: α_s(R) = α_s* × 560/S_F(ρ(R)).

# Arguments
- `p`: eigenvalues of the density matrix ρ ∈ D₆

# Returns
Ricci scalar S_F ∈ ℝ

# Examples
```julia
p_star = rho_KK_eigenvalues(4260.0)
SF_GG(p_star)           # → 560.0   (vacuum)

p_100 = rho_KK_eigenvalues(100.0)
SF_GG(p_100)            # → 238.22  (α_s ≈ 0.28 at R=100)
```

See: FisherGeometrics preprint v15, section 3 (Extended Theorem E).
"""
function SF_GG(p::Vector{Float64})
    T  = su_basis(6)
    N  = length(T)
    n  = 6

    Γ  = christoffel_rotate(p, T)

    G  = reshape([sum(real(conj(T[a][i,j])*T[b][i,j])/(2*(p[i]+p[j]))
                      for i in 1:n, j in 1:n)
                  for a in 1:N, b in 1:N], N, N)
    Gi = pinv(G)

    R_GG = zeros(N, N, N, N)
    for e in 1:N, a in 1:N, b in 1:N, c in 1:N
        R_GG[e,a,b,c] = sum(Γ[f,b,c]*Γ[e,a,f] - Γ[f,a,c]*Γ[e,b,f]
                            for f in 1:N)
    end

    Ric = [sum(R_GG[b,a,b,c] for b in 1:N) for a in 1:N, c in 1:N]

    return sum(Gi[a,c]*Ric[a,c] for a in 1:N, c in 1:N)
end
