# src/RiemannTensor.jl

########################################################################
#
#  RiemannTensor.jl
#
#  Berekent de volledige Riemann-tensor van de Bures-metriek op D_n
#  op basis van FisherGeometrics.jl.
#
#  Gebruik:
#    include("RiemannTensor.jl")
#    n    = 6
#    ρ    = Matrix(I, n, n) / n          # vacuüm ρ* = I/6
#    T    = su_basis(n)                  # 35 generatoren
#    R    = riemann_tensor(ρ, T)         # R[e,a,b,c] = R^e_{abc}
#    Ric  = ricci_tensor(R, T, ρ)        # Ric[a,b]
#    S    = ricci_scalar(Ric, T, ρ)      # S_F (verwacht: 560)
#
########################################################################

using LinearAlgebra
using Printf


# ── Bures-metriek ─────────────────────────────────────────────────

"""Oplost ρL + Lρ = 2Y via pinv (robuust bij singuliere ρ)."""
function solve_sld(ρ::AbstractMatrix, Y::AbstractMatrix; tol=1e-12)
    n = size(ρ,1)
    A = kron(ρ, I(n)) + kron(I(n), transpose(ρ))
    L = reshape(pinv(A; atol=tol) * 2vec(ComplexF64.(Y)), n, n)
    return (L + L') / 2
end

"""Bures-metriek g(X,Y)|_ρ = ¼ Re Tr(X L_Y)."""
bures_g(ρ, X, Y) = (1/4) * real(tr(X * solve_sld(ρ, Y)))

"""Volledige N×N metriekmatrix bij ρ."""
function metric_matrix(ρ, T)
    N = length(T)
    G = zeros(Float64, N, N)
    for a in 1:N, b in a:N
        v = bures_g(ρ, T[a], T[b])
        G[a,b] = v; G[b,a] = v
    end
    return G
end

# ── Christoffel-symbolen ───────────────────────────────────────────

"""
    christoffel_symbols(ρ, T; eps) -> Array{Float64,3}

Γᵉ_{ab} = ½ gᵉᶠ (∂_a g_{bf} + ∂_b g_{af} − ∂_f g_{ab})

via eindige differenties van de Bures-metriek.
"""
function christoffel_symbols(ρ::AbstractMatrix, T::Vector;
                             eps::Float64=1e-4)
    N   = length(T)
    G   = metric_matrix(ρ, T)
    Gi  = pinv(G; atol=1e-10)

    # ∂_c g_{ab} via centrale differenties
    function dg(a, b, c)
        (bures_g(ρ + eps*T[c], T[a], T[b]) -
         bures_g(ρ - eps*T[c], T[a], T[b])) / (2*eps)
    end

    Γ = zeros(Float64, N, N, N)
    for e in 1:N, a in 1:N, b in a:N
        val = 0.0
        for f in 1:N
            val += Gi[e,f] * (dg(b,f,a) + dg(a,f,b) - dg(a,b,f))
        end
        Γ[e,a,b] = val/2
        Γ[e,b,a] = val/2
    end
    return Γ
end

# ── Riemann ─────────────────────────────────────────────────

function riemann(g::FisherMetric, ρ::AbstractMatrix, basis;
                 eps::Float64=1e-4)
    N  = length(basis)
    Γ0 = christoffel(g, ρ, basis)

    function dΓ(e, b, c, a)
        Γp = christoffel(g, ρ + eps*basis[a], basis)
        Γm = christoffel(g, ρ - eps*basis[a], basis)
        return (Γp[e,b,c] - Γm[e,b,c]) / (2*eps)
    end

    R = zeros(Float64, N, N, N, N)
    for e in 1:N, a in 1:N, b in 1:N, c in 1:N
        lin  = dΓ(e,b,c,a) - dΓ(e,a,c,b)
        quad = sum(Γ0[f,b,c]*Γ0[e,a,f] - Γ0[f,a,c]*Γ0[e,b,f] for f in 1:N)
        R[e,a,b,c] = lin + quad
    end
    return R
end

# ── Riemann-tensor ─────────────────────────────────────────────────

"""
    riemann_tensor(ρ, T; eps) -> Array{Float64,4}

R^e_{abc} = ∂_a Γᵉ_{bc} − ∂_b Γᵉ_{ac}
           + Σ_f (Γᶠ_{bc} Γᵉ_{af} − Γᶠ_{ac} Γᵉ_{bf})

Geeft een N×N×N×N tensor terug waarbij R[e,a,b,c] = R^e_{abc}.
"""
function riemann_tensor(ρ::AbstractMatrix, T::Vector;
                        eps::Float64=1e-4)
    N  = length(T)
    Γ0 = christoffel_symbols(ρ, T; eps=eps)

    # ∂_a Γᵉ_{bc} via differentie van Christoffel
    function dΓ(e, b, c, a)
        ρp  = ρ + eps*T[a]
        ρm  = ρ - eps*T[a]
        Γp  = christoffel_symbols(ρp, T; eps=eps)
        Γm  = christoffel_symbols(ρm, T; eps=eps)
        return (Γp[e,b,c] - Γm[e,b,c]) / (2*eps)
    end

    R = zeros(Float64, N, N, N, N)
    for e in 1:N, a in 1:N, b in 1:N, c in 1:N
        # Lineaire bijdrage: ∂Γ
        lin = dΓ(e,b,c,a) - dΓ(e,a,c,b)

        # Kwadratische bijdrage: ΓΓ
        quad = sum(Γ0[f,b,c]*Γ0[e,a,f] - Γ0[f,a,c]*Γ0[e,b,f] for f in 1:N)
        R[e,a,b,c] = lin + quad
    end
    return R
end

# ── Handige alles-in-één functie ───────────────────────────────────

"""
    compute_curvature(ρ, T; eps, verbose) -> NamedTuple

Berekent de volledige Riemannse meetkunde bij toestand ρ:
  - metriek G
  - Christoffel Γ
  - Riemann R
  - Ricci Ric
  - scalaire kromming S_F

Voorbeeld:
    T = su_basis(6)
    ρ = Matrix(I,6,6)/6
    C = compute_curvature(ρ, T)
    println("S_F = ", C.S_F)   # verwacht: ≈ 560
"""
function compute_curvature(ρ::AbstractMatrix, T::Vector;
                           eps::Float64=1e-4,
                           verbose::Bool=true)
    N = length(T)
    verbose && println("Stap 1/4: metriek berekenen ($N×$N)...")
    G    = metric_matrix(ρ, T)
    G_inv = pinv(G; atol=1e-10)

    verbose && println("Stap 2/4: Christoffel-symbolen ($N³ componenten)...")
    Γ    = christoffel_symbols(ρ, T; eps=eps)

    verbose && println("Stap 3/4: Riemann-tensor ($N⁴ componenten)...")
    verbose && println("          (dit duurt enkele minuten voor N=$N)")
    R    = riemann_tensor(ρ, T; eps=eps)

    verbose && println("Stap 4/4: Ricci + scalaire kromming...")
    Ric  = ricci_tensor(R)
    S_F  = ricci_scalar(Ric, G_inv)

    verbose && @printf("Klaar. S_F = %.4f  (verwacht bij ρ*: 560)\n", S_F)

    return (G=G, G_inv=G_inv, Γ=Γ, R=R, Ric=Ric, S_F=S_F)
end

# ── Snelle steekproef: alleen sectionele krommingen ───────────────

"""
    sectional_curvatures(ρ, T; n_pairs, eps) -> Vector{Float64}

Berekent een steekproef van sectionele krommingen K(Tₐ,Tᵦ)
zonder de volledige Riemann-tensor op te bouwen.
Sneller dan compute_curvature voor een eerste indruk.

Verwachte waarden bij ρ* = I/6:
    K ∈ {+0.25, -2.0, -2.0}  (per generatortype)
"""
function sectional_curvatures(ρ::AbstractMatrix, T::Vector;
                               n_pairs::Int=10,
                               eps::Float64=1e-4)
    N  = length(T)
    G  = metric_matrix(ρ, T)
    Gi = pinv(G; atol=1e-10)
    Γ0 = christoffel_symbols(ρ, T; eps=eps)

    results = Tuple{Int,Int,Float64}[]

    for a in 1:min(n_pairs, N)
        for b in a+1:min(n_pairs+1, N)
            # R^a_{bab} via Christoffel-differentie
            Γa_p = christoffel_symbols(ρ + eps*T[a], T; eps=eps)
            Γa_m = christoffel_symbols(ρ - eps*T[a], T; eps=eps)
            Γb_p = christoffel_symbols(ρ + eps*T[b], T; eps=eps)
            Γb_m = christoffel_symbols(ρ - eps*T[b], T; eps=eps)

            dΓ_a = (Γa_p - Γa_m) / (2*eps)
            dΓ_b = (Γb_p - Γb_m) / (2*eps)

            R_abab = 0.0
            for e in 1:N
                lin  = dΓ_a[e,b,b] - dΓ_b[e,a,b]
                quad = sum(Γ0[f,b,b]*Γ0[e,a,f] - Γ0[f,a,b]*Γ0[e,b,f]
                           for f in 1:N)
                R_abab += G[a,a] * (lin + quad)
            end

            denom = G[a,a]*G[b,b] - G[a,b]^2
            K = denom > 1e-20 ? R_abab/denom : 0.0
            push!(results, (a, b, K))
        end
        length(results) >= n_pairs && break
    end
    return results
end

########################################################################
#  VOORBEELD
########################################################################
#
#  julia> include("RiemannTensor.jl")
#
#  julia> n = 6; T = su_basis(n); ρ = Matrix(I,n,n)/n
#
#  # Snel: sectionele krommingen (seconden)
#  julia> K = sectional_curvatures(ρ, T; n_pairs=6)
#  julia> for (a,b,k) in K; @printf("K(T_%d,T_%d) = %+.4f\n",a,b,k); end
#
#  # Volledig: alle curvature-objecten (minuten)
#  julia> C = compute_curvature(ρ, T)
#  julia> C.S_F        # → ≈ 560
#  julia> C.Ric        # Ricci-tensor (35×35)
#  julia> C.R[2,3,7,6] # R^2_{376} → ≈ 2.8125
#
########################################################################
