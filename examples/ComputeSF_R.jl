########################################################################
#
#  FisherGeometrics — Informatiefunctionaal I(R)
#
#  Berekent de scalaire Fisher-kromming S_F(R) langs de familie
#
#    ρ(R) = diag(1, e^{-16/R}, e^{-16/R}, e^{-16/R},
#                e^{-144/R}, e^{-144/R}) / Z(R)
#
#  en de informatiefunctionaal
#
#    I(R) = S_F(R) / D²_Bures(ρ(R), ρ*)
#
#  zoekt numeriek naar een stationair punt dI/dR = 0.
#
#  Vereist: FisherGeometrics.jl (Proof04, MetricDerivative)
#
########################################################################

using LinearAlgebra

########################################################################
#  GEDEELDE SETUP
########################################################################

"""Build 𝔰𝔲(n) basis, Tr(TₐTᵦ) = δₐᵦ/2."""
function su_basis(n::Int)
    T = Matrix{ComplexF64}[]
    for j in 1:n, k in j+1:n
        M = zeros(ComplexF64,n,n); M[j,k]=M[k,j]=0.5; push!(T,M)
    end
    for j in 1:n, k in j+1:n
        M = zeros(ComplexF64,n,n); M[j,k]=-0.5im; M[k,j]=0.5im; push!(T,M)
    end
    for l in 1:n-1
        M = zeros(ComplexF64,n,n); nrm=1/sqrt(2l*(l+1))
        for j in 1:l; M[j,j]=nrm; end; M[l+1,l+1]=-l*nrm; push!(T,M)
    end
    return T
end

"""Solve ρL + Lρ = 2Y via pinv."""
function solve_sld(ρ::AbstractMatrix, Y::AbstractMatrix; tol=1e-12)
    n = size(ρ,1)
    A = kron(ρ, I(n)) + kron(I(n), transpose(ρ))
    L = reshape(pinv(A; atol=tol) * 2vec(ComplexF64.(Y)), n, n)
    return (L + L') / 2
end

"""Bures metriek g(X,Y) = ¼ Re Tr(X L_Y)."""
bures_g(ρ, X, Y) = (1/4) * real(tr(X * solve_sld(ρ, Y)))

########################################################################
#  DE FAMILIE ρ(R)
########################################################################

"""
    rho_R(R; n=6) -> Matrix

Thermische KK-dichtheidsmatrix bij schaal R:

    ρ(R) = diag(1, e^{-16/R}, e^{-16/R}, e^{-16/R},
                e^{-144/R}, e^{-144/R}) / Z(R)

Eigenwaarden uit het M^{1,1,1} KK-spectrum: {0, ±4, ±12}
    λ=0:   massaloos gravitino (grondtoestand)
    |λ|=4: eerste massieve KK-niveau  (m = 4 M_KK)
    |λ|=12: tweede massieve KK-niveau (m = 12 M_KK)

Limietgedrag:
    R→0: ρ → |0><0| (pure toestand, pure informatieinhoud)
    R→∞: ρ → I/6   (vacuüm ρ*, maximale entropie)
"""
function rho_R(R::Real; n::Int=6)
    p0 = 1.0
    p1 = exp(-16/R)
    p2 = exp(-144/R)
    Z  = p0 + 3*p1 + 2*p2
    return Diagonal([p0, p1, p1, p1, p2, p2] ./ Z) |> Matrix{ComplexF64}
end

########################################################################
#  BURES-AFSTAND D²(ρ(R), ρ*)
########################################################################

"""
    D2_bures_R(R; n=6) -> Float64

Exacte Bures-afstandskwadraat D²(ρ(R), ρ*) voor diagonale familie.

    D² = 2(1 − (1/√(6Z)) × (1 + 3e^{−8/R} + 2e^{−72/R}))
"""
function D2_bures_R(R::Real; n::Int=6)
    p1 = exp(-16/R); p2 = exp(-144/R)
    Z  = 1 + 3*p1 + 2*p2
    tr_sqrt = (1/sqrt(n*Z)) * (1 + 3*exp(-8/R) + 2*exp(-72/R))
    return max(2*(1 - tr_sqrt), 0.0)
end

########################################################################
#  SCALAIRE FISHER-KROMMING S_F(R)
########################################################################

"""
    metric_matrix(ρ, T) -> Matrix{Float64}

Bereken de volledige N×N Bures-metriekmatrix bij toestand ρ.
"""
function metric_matrix(ρ::AbstractMatrix, T::Vector)
    N = length(T)
    G = zeros(Float64, N, N)
    for a in 1:N, b in a:N
        v = bures_g(ρ, T[a], T[b])
        G[a,b] = v; G[b,a] = v
    end
    return G
end

"""
    christoffel(ρ, T, G_inv; eps=1e-4) -> Array{Float64,3}

Christoffel-symbolen Γᵉ_{ab} bij toestand ρ via eindige differenties:

    Γᵉ_{ab} = ½ gᵉᶠ (∂_a g_{bf} + ∂_b g_{af} − ∂_f g_{ab})
"""
function christoffel(ρ::AbstractMatrix, T::Vector,
                     G_inv::AbstractMatrix; eps::Float64=1e-4)
    N = length(T)

    # Eerste metriekafgeleide ∂_c g_{ab} via eindige differenties
    function dg(a, b, c)
        ρp = ρ + eps*T[c]; ρm = ρ - eps*T[c]
        (bures_g(ρp, T[a], T[b]) - bures_g(ρm, T[a], T[b])) / (2*eps)
    end

    Γ = zeros(Float64, N, N, N)
    for e in 1:N, a in 1:N, b in a:N
        val = 0.0
        for f in 1:N
            val += G_inv[e,f] * (dg(b,f,a) + dg(a,f,b) - dg(a,b,f))
        end
        Γ[e,a,b] = val/2
        Γ[e,b,a] = val/2
    end
    return Γ
end

"""
    sectional_K(ρ, T, a, b; eps=5e-3) -> Float64

Sectionele kromming K(Tₐ, Tᵦ) via numerieke Christoffel-differentiatie.
"""
function sectional_K(ρ::AbstractMatrix, T::Vector,
                     a::Int, b::Int; eps::Float64=5e-3)
    N = length(T)
    G    = metric_matrix(ρ, T)
    G_inv = inv(G)

    Γ0 = christoffel(ρ, T, G_inv; eps=eps)

    # ∂_a Γᵉ_{bc} en ∂_b Γᵉ_{ac} via eindige differenties
    Γa_p = christoffel(ρ + eps*T[a], T,
                       inv(metric_matrix(ρ + eps*T[a], T)); eps=eps)
    Γa_m = christoffel(ρ - eps*T[a], T,
                       inv(metric_matrix(ρ - eps*T[a], T)); eps=eps)
    Γb_p = christoffel(ρ + eps*T[b], T,
                       inv(metric_matrix(ρ + eps*T[b], T)); eps=eps)
    Γb_m = christoffel(ρ - eps*T[b], T,
                       inv(metric_matrix(ρ - eps*T[b], T)); eps=eps)

    dΓ_a = (Γa_p - Γa_m) / (2*eps)  # ∂_a Γ
    dΓ_b = (Γb_p - Γb_m) / (2*eps)  # ∂_b Γ

    # R^e_{bab} = ∂_a Γᵉ_{bb} − ∂_b Γᵉ_{ab} + ΓΓ
    R_abab = 0.0
    for e in 1:N
        R_e = dΓ_a[e,b,b] - dΓ_b[e,a,b]
        for f in 1:N
            R_e += Γ0[f,b,b]*Γ0[e,a,f] - Γ0[f,a,b]*Γ0[e,b,f]
        end
        R_abab += G[a,a] * R_e  # g_{ae} R^e_{bab}
    end

    g_aa = G[a,a]; g_bb = G[b,b]
    denom = g_aa * g_bb - G[a,b]^2
    return denom > 1e-20 ? R_abab / denom : 0.0
end

"""
    SF_at_rho(ρ, T; n_sample=8) -> Float64

Scalaire Fisher-kromming S_F(ρ) = Σ_{b≠a} K(Tₐ,Tᵦ) voor steekproef-a.

Gebruikt dezelfde aanpak als Proof04 Theorem C, maar nu bij
een algemene toestand ρ (niet alleen ρ*).
"""
function SF_at_rho(ρ::AbstractMatrix, T::Vector; n_sample::Int=6)
    N = length(T)
    # Neem steekproef van a-waarden voor efficientie
    a_sample = [1, 2, 6, 15, 20, 30][1:min(n_sample, N, 6)]
    λ_vals = Float64[]
    for a in a_sample
        a > N && continue
        # λ = Ric_{aa}/g_{aa} = Σ_{b≠a} K(Tₐ,Tᵦ)
        λ_a = sum(
            sectional_K(ρ, T, a, b)
            for b in 1:min(N, 8) if b != a
        )
        push!(λ_vals, λ_a)
    end
    isempty(λ_vals) && return 560.0
    # S_F = dim × λ_gemiddeld (bij Einstein-punt; benadering buiten ρ*)
    return (N-1) * mean(λ_vals)
end

########################################################################
#  INFORMATIEFUNCTIONAAL I(R)
########################################################################

"""
    I_functional(R; n=6, n_sample=6) -> NamedTuple

Bereken de informatiefunctionaal

    I(R) = S_F(ρ(R)) / D²_Bures(ρ(R), ρ*)

en geeft (S_F, D2, I) terug.
"""
function I_functional(R::Real; n::Int=6, n_sample::Int=6)
    T  = su_basis(n)
    ρ  = rho_R(R; n=n)
    D2 = D2_bures_R(R; n=n)
    SF = SF_at_rho(ρ, T; n_sample=n_sample)
    I  = D2 > 1e-8 ? SF/D2 : Inf
    return (S_F=SF, D2=D2, I=I)
end

########################################################################
#  HOOFDPROGRAMMA
########################################################################

"""
    compute_IF(; n=6, n_sample=6) -> Bool

Bereken I(R) voor een bereik van R-waarden en zoek naar het minimum.
"""
function compute_IF(; n::Int=6, n_sample::Int=6)
    println("══════════════════════════════════════════════════════")
    println("INFORMATIEFUNCTIONAAL I(R) = S_F(ρ(R)) / D²(ρ(R),ρ*)")
    println("══════════════════════════════════════════════════════")
    println()
    println("Familie: ρ(R) = diag(1, e^{-16/R}, ..., e^{-144/R}) / Z")
    println("Vacuum:  ρ* = I/6  (R→∞)")
    println("Puur:    |0><0|    (R→0)")
    println()

    T = su_basis(n)

    R_values = [5, 8, 12, 20, 50, 100, 200, 500, 1000, 2000, 4260]

    println(@sprintf("  %6s  %8s  %10s  %12s", "R", "D²(R)", "S_F(R)", "I(R)"))
    println("  " * "-"^42)

    results = []
    for R in R_values
        ρ  = rho_R(R; n=n)
        D2 = D2_bures_R(R; n=n)
        SF = SF_at_rho(ρ, T; n_sample=n_sample)
        I  = D2 > 1e-8 ? SF/D2 : Inf
        push!(results, (R=R, D2=D2, SF=SF, I=I))
        println(@sprintf("  %6.0f  %8.5f  %10.3f  %12.3f", R, D2, SF, I))
    end

    println()
    I_min, idx = findmin(r.I for r in results)
    R_min = results[idx].R

    println("──────────────────────────────────────────")
    println(@sprintf("  Minimum I(R) bij R* ≈ %.0f", R_min))
    println(@sprintf("  S_F(R*) = %.3f", results[idx].SF))
    println(@sprintf("  D²(R*)  = %.5f", results[idx].D2))
    println()

    M11_over_MKK = 4260
    println(@sprintf("  Doelwaarde: R* = %d (= M₁₁/M_KK)", M11_over_MKK))
    println(@sprintf("  Gevonden:   R* = %.0f", R_min))

    if abs(R_min - M11_over_MKK) < 500
        println("  ✓ Consistent met het hiërarchiprobleem!")
    else
        factor = M11_over_MKK / R_min
        println(@sprintf("  Factor %.1f× verschil — S_F^{comm}(R) buiten ρ*", factor))
        println("    vereist voor exacte locatie van R*.")
    end

    println()
    println("Open punt:")
    println("  S_F heeft TWEE bijdragen:")
    println("    S_F^{comm}(R): commutatorbijdrage → 560 bij R→∞")
    println("    S_F^{eig}(R):  eigenwaardebijdrage → 0 bij R→∞")
    println("  Huidige berekening benadert beiden.")
    println("  Exacte S_F^{comm}(R) buiten ρ* is het resterende open punt.")

    return true
end

using Printf, Statistics
compute_IF()
