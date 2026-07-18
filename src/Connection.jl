# src/Connection.jl
 
abstract type AbstractConnection end
 
struct LeviCivitaConnection <: AbstractConnection
end
 
"""
    christoffel(g, ρ, basis) -> Array{Float64,3}
 
Levi-Civita Christoffel-symbolen van de Fisher/Bures-metriek.
 
Conventie:  Γ[e,a,b] = Γᵉₐᵦ  (bovenste index EERSTE)
 
    ∇_{eₐ} eᵦ = Σₑ Γᵉₐᵦ eₑ
 
    Γᵉₐᵦ = ½ gᵉˡ (∂ₐ g_{bℓ} + ∂ᵦ g_{aℓ} − ∂ℓ g_{ab})
 
Consistent met `christoffel_vacuum` in Curvature.jl en
`christoffel_rotate` hieronder.
 
!!! note "Index-conventie"
    Vóór juli 2026 gebruikte deze functie Γ[i,j,k] = Γᵏᵢⱼ
    (bovenste index derde). De conventie is omgekeerd om consistent
    te zijn met de rest van de codebase.
"""
function christoffel(g::FisherMetric, ρ, basis)
    G    = metric_matrix(g, ρ, basis)
    Ginv = pinv(G; atol=1e-10)
    dg   = metric_derivatives(g, ρ, basis)
    N    = length(basis)
 
    Γ = zeros(eltype(G), N, N, N)
    for e in 1:N, a in 1:N, b in 1:N
        Γ[e,a,b] = 0.5 * sum(
            Ginv[e,l] * (dg[b,l,a] + dg[a,l,b] - dg[a,b,l])
            for l in 1:N
        )
    end
    return Γ
end
 
"""
    christoffel_rotate(p, basis) -> Array{Float64,3}
 
Exacte Christoffel-symbolen voor diagonale ρ = diag(p) via
`dmetric_rotate`.
 
Conventie:  Γ[e,a,b] = Γᵉₐᵦ  (bovenste index EERSTE)
 
Voordelen t.o.v. `christoffel()`:
- Exact: geen finite differences van de metriek
- Correct voor alle R, ook kleine R en degenerate eigenwaarden (pₗ = pₘ)
- ~10× sneller door directe SLD-berekening
 
# Argumenten
- `p`: eigenwaarden van ρ als `Vector{Float64}`
- `basis`: vector van su(n)-generatoren
 
# Gebruik
```julia
p = rho_KK_eigenvalues(100.0)
T = su_basis(6)
Γ = christoffel_rotate(p, T)   # 35×35×35 tensor, Γ[e,a,b] = Γᵉₐᵦ
```
 
# Vergelijk
- `christoffel()`: via `metric_derivatives` (finite diff van bures_g)
- `christoffel_vacuum()`: exact analytisch bij ρ* = I/n via d-symbolen
"""
function christoffel_rotate(p::Vector{<:Real}, basis)
    N  = length(basis)
    n  = length(p)
 
    # Metriek
    G  = zeros(Float64, N, N)
    for a in 1:N, b in 1:N
        G[a,b] = sum(real(conj(basis[a][i,j])*basis[b][i,j])/(2*(p[i]+p[j]))
                     for i in 1:n, j in 1:n)
    end
    Gi = pinv(G; atol=1e-10)
 
    # Exacte metriekafgeleiden
    dg = metric_derivatives_rotate(p, basis)
 
    # Christoffel — zelfde formule als christoffel(), zelfde conventie
    Γ = zeros(Float64, N, N, N)
    for e in 1:N, a in 1:N, b in 1:N
        Γ[e,a,b] = 0.5 * sum(
            Gi[e,l] * (dg[b,l,a] + dg[a,l,b] - dg[a,b,l])
            for l in 1:N
        )
    end
    return Γ
end
