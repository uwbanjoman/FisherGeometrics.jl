# src/MetricDerivative.jl

using LinearAlgebra

"""
    dmetric(metric, ρ, X, Y, Z)

Directional derivative of the Fisher/Bures metric
along the tangent direction `Z`.

Computes

    D_Z gρ(X,Y)

using

    D(Lρ⁻¹)
        =
        -Lρ⁻¹ ∘ D(Lρ) ∘ Lρ⁻¹.
"""
function dmetric(::FisherMetric,
                 ρ::AbstractMatrix,
                 X::AbstractMatrix,
                 Y::AbstractMatrix,
                 Z::AbstractMatrix)

    # Lρ⁻¹(Y)
    LY = Lρ_inv(ρ, Y)

    # D(Lρ)(LY)
    DLY = dLρ(Z, LY)

    # Lρ⁻¹(D(Lρ)(LY))
    term = Lρ_inv(ρ, DLY)

    # minus sign from
    #
    # D(L⁻¹)=−L⁻¹(DL)L⁻¹
    #
    #return -real(tr(X * term))/2
    return -real(tr(X * term))/4

end

"""
    metric_derivatives(g, ρ, basis)

Compute the directional derivatives of the Fisher metric in the
basis `basis`.

Returns a 3-tensor `dg` with

    dg[i,j,k] = Dgρ(H_k; X_i, X_j),

where

- `X_i = basis[i]`,
- `H_k = basis[k]`.

The tensor `dg` is symmetric in its first two indices.
"""
function metric_derivatives(g::FisherMetric,
                            ρ::AbstractMatrix,
                            basis)
    n = length(basis)

    dg = zeros(Float64, n, n, n)

    for i in 1:n
        Xi = basis[i]

        for j in i:n
            Xj = basis[j]

            for k in 1:n
                Hk = basis[k]

                val = dmetric(g, ρ, Xi, Xj, Hk)

                dg[i,j,k] = val
                dg[j,i,k] = val
            end
        end
    end

    return dg
end

"""
    metric_derivatives(g::FisherMetric, ρ::Diagonal, basis)

Geoptimaliseerde directionele afgeleiden voor diagonale density matrices.
Voorkomt de zware overhead van numerieke differentiatie.
"""
function metric_derivatives(g::FisherMetric, ρ::Diagonal{T}, basis) where {T<:Real}
    n = length(basis)
    dg = zeros(Float64, n, n, n)
    mat_dim = size(ρ, 1)
    diag_ρ = ρ.diag

    # Pre-allocateer de kwadratische denominators
    inv_denom_sq = Matrix{T}(undef, mat_dim, mat_dim)
    for j in 1:mat_dim
        λ_j = diag_ρ[j]
        for i in 1:mat_dim
            denom = diag_ρ[i] + λ_j
            inv_denom_sq[i, j] = denom > 1e-12 ? 1.0 / (denom^2) : 0.0
        end
    end

    # Als we op een maximally mixed state zitten, zijn de afgeleiden 0.0.
    # We berekenen het hier expliciet via de analytische element-wise formule:
    for i in 1:n
        Xi = basis[i]
        for j in i:n
            Xj = basis[j]
            for k in 1:n
                Hk = basis[k]
                
                accum = 0.0
                @inbounds for b in 1:mat_dim
                    for a in 1:mat_dim
                        # Richtingafgeleide-bijdrage van de diagonale verschuiving
                        # (Hk[a,a] + Hk[b,b]) vertegenwoordigt dρ
                        dH = real(Hk[a, a] + Hk[b, b])
                        accum += real(Xi[b, a] * Xj[a, b]) * dH * inv_denom_sq[a, b]
                    end
                end
                
                val = -0.5 * accum
                dg[i, j, k] = val
                dg[j, i, k] = val
            end
        end
    end

    return dg
end

"""
    ddmetric(g, ρ, X, Y, H, K)

Second directional derivative of the Fisher metric.

Returns

    D²gρ(H,K)(X,Y)

where `H` and `K` are tangent directions and `X`,`Y` are the metric
arguments.

The implementation uses

    D²Lρ = 0

since the Jordan operator is linear in `ρ`.
"""
function ddmetric(::FisherMetric,
                  ρ::AbstractMatrix,
                  X::AbstractMatrix,
                  Y::AbstractMatrix,
                  H::AbstractMatrix,
                  K::AbstractMatrix)

    # Lρ⁻¹(Y)
    LY = Lρ_inv(ρ, Y)

    #
    # First contribution
    #
    T1 = Lρ_inv(ρ,
                jordan(K,
                    Lρ_inv(ρ,
                        jordan(H, LY))))

    #
    # Second contribution
    #
    T2 = Lρ_inv(ρ,
                jordan(H,
                    Lρ_inv(ρ,
                        jordan(K, LY))))

    #return real(tr(X * (T1 + T2))) / 2
    return real(tr(X * (T1 + T2)))/4

end

"""
    ddmetric_tensor(g, ρ, basis)

Compute the second directional derivatives of the Fisher metric in the
basis of tangent vectors.

Returns the rank-4 tensor

    ddg[i,j,k,l] =
        D²gρ(Hₖ,Hₗ)(Xᵢ,Xⱼ),

where

- `Xᵢ = basis[i]`,
- `Hₖ = basis[k]`.

The tensor is symmetric in

- `(i,j)` because the metric is symmetric;
- `(k,l)` because mixed directional derivatives commute.
"""
function ddmetric_tensor(g::FisherMetric,
                         ρ::AbstractMatrix,
                         basis)
    n = length(basis)
    T = zeros(Float64, n, n, n, n)

    for i in 1:n
        Xi = basis[i]

        for j in 1:n
            Yj = basis[j]

            for k in 1:n
                Hk = basis[k]

                for l in 1:n
                    Hl = basis[l]

                    T[i,j,k,l] =
                        ddmetric(g,
                                 ρ,
                                 Xi,
                                 Yj,
                                 Hk,
                                 Hl)

                end
            end
        end
    end

    return T

end

"""
    ddmetric_tensor(g::FisherMetric, ρ::Diagonal, basis)

Geoptimaliseerde tweede directionele afgeleiden voor diagonale density matrices.
Elimineert de honderdduizenden numerieke differentiatie-allocaties.
"""
function ddmetric_tensor(g::FisherMetric, ρ::Diagonal{T}, basis) where {T<:Real}
    n = length(basis)
    ddg = zeros(Float64, n, n, n, n)
    mat_dim = size(ρ, 1)
    diag_ρ = ρ.diag

    # Pre-allocateer de derdegraads denominators (λ_i + λ_j)³
    inv_denom_cube = Matrix{T}(undef, mat_dim, mat_dim)
    for j in 1:mat_dim
        λ_j = diag_ρ[j]
        for i in 1:mat_dim
            denom = diag_ρ[i] + λ_j
            inv_denom_cube[i, j] = denom > 1e-12 ? 1.0 / (denom^3) : 0.0
        end
    end

    # Analytische 4D contractielus
    for i in 1:n
        Xi = basis[i]
        for j in i:n
            Xj = basis[j]
            
            for k in 1:n
                Hk = basis[k]
                for l in k:n
                    Hl = basis[l]
                    
                    accum = 0.0
                    @inbounds for b in 1:mat_dim
                        for a in 1:mat_dim
                            # Tweede orde afgeleide-bijdrage van de diagonale verschuivingen
                            dH_k = real(Hk[a, a] + Hk[b, b])
                            dH_l = real(Hl[a, a] + Hl[b, b])
                            
                            accum += real(Xi[b, a] * Xj[a, b]) * dH_k * dH_l * inv_denom_cube[a, b]
                        end
                    end
                    
                    # De tweede afgeleide van 1/(p_i+p_j) levert een factor 2 * (1/(p_i+p_j)³) 
                    # Gecombineerd met de 0.5 van de metriek-definitie geeft dit val:
                    val = accum # (0.5 * 2 = 1.0)
                    
                    # Profiteer van alle symmetrieën in de 4-tensor:
                    # Symmetrisch in (i,j) en symmetrisch in (k,l)
                    ddg[i, j, k, l] = val
                    ddg[j, i, k, l] = val
                    ddg[i, j, l, k] = val
                    ddg[j, i, l, k] = val
                end
            end
        end
    end

    return ddg
end

# ── Exacte metriekafgeleide voor diagonale toestanden ─────────────
 
"""
    dmetric_rotate(p, X, Y, Z) -> Float64
 
Exacte richtingsafgeleide van de Bures-metriek bij diagonale ρ = diag(p):
 
    D_Z g(X,Y)|_ρ
 
via de SLD-formule (alle generatoren) gecombineerd met de
eigenvector-rotatie-interpretatie (Haber 2019) voor off-diagonale Z.
 
Drie gevallen op basis van het type Z:
 
1. **Diagonale Z** (Z[i,j]=0 voor i≠j):
   Eigenwaarden van ρ veranderen, eigenvectoren niet.
   Exact analytische formule via (c_i+c_j)/(p_i+p_j)².
 
2. **Off-diagonale Z** (symmetrisch of antisymmetrisch):
   Eigenvectoren roteren bij perturbatie ρ+εZ.
   SLD-formule: D_Z g(X,Y) = −½ Re Tr(X × Lρ⁻¹(jordan(Z, Lρ⁻¹(Y))))
   met Lρ⁻¹(Y)[i,j] = Y[i,j]/(p_i+p_j).
 
Geldig voor **alle R** — ook kleine R met grote eigenwaardenverschillen,
en ook degenerate eigenwaarden (p_l = p_m).
 
# Argumenten
- `p`: eigenwaarden van ρ als Vector{Float64}
- `X`, `Y`: Hermitische generatoren (metriek-argumenten)
- `Z`: Hermitische generator (differentiatierichting)
 
# Vergelijk met
- `dmetric(::FisherMetric, ρ, X, Y, Z)`: via `Lρ_inv` (volle matrix)
- `metric_derivatives(::FisherMetric, ρ, basis)`: volledige tensor
 
# Gebruik
```julia
p   = rho_KK_eigenvalues(100.0)
T   = su_basis(6)
val = dmetric_rotate(p, T[1], T[1], T[32])   # ≈ -1.299
```
"""
function dmetric_rotate(p::Vector{<:Real},
                        X::AbstractMatrix,
                        Y::AbstractMatrix,
                        Z::AbstractMatrix)
    n = length(p)
 
    # Bepaal type Z
    is_diag = all(abs(Z[i,j]) < 1e-12 for i in 1:n, j in 1:n if i ≠ j)
 
    if is_diag
        # Geval 1: diagonale Z — eigenwaarden veranderen
        # D_Z g(X,Y) = −½ Σ_{ij} Re(X*_{ij} Y_{ij}) (c_i+c_j)/(p_i+p_j)²
        s = 0.0
        for i in 1:n, j in 1:n
            d = p[i]+p[j]; d < 1e-15 && continue
            ci = real(Z[i,i]); cj = real(Z[j,j])
            s -= real(conj(X[i,j])*Y[i,j]) * (ci+cj) / d^2
        end
        return s / 2
 
    else
        # Geval 2/3: off-diagonale Z — eigenvector-rotatie
        # D_Z g(X,Y) = −½ Re Tr(X × Lρ⁻¹(jordan(Z, Lρ⁻¹(Y))))
        # met Lρ⁻¹(A)[i,j] = A[i,j]/(p_i+p_j)
        LY   = [Y[i,j]/(p[i]+p[j]) for i in 1:n, j in 1:n]
        jord = Z * LY + LY * Z
        LjY  = [jord[i,j]/(p[i]+p[j]) for i in 1:n, j in 1:n]
        return -real(tr(X * LjY)) / 2
    end
end
 
"""
    metric_derivatives_rotate(p, basis) -> Array{Float64,3}
 
Volledige tensor van metriekafgeleiden via `dmetric_rotate`:
 
    dg[a,b,c] = D_{basis[c]} g(basis[a], basis[b])|_{diag(p)}
 
Exact voor alle generatoren, ook off-diagonale en degenerate eigenwaarden.
Sneller dan `metric_derivatives` via SLD-inverse voor diagonale ρ.
 
# Gebruik
```julia
p  = rho_KK_eigenvalues(100.0)
T  = su_basis(6)
dg = metric_derivatives_rotate(p, T)   # 35×35×35 tensor
```
"""
function metric_derivatives_rotate(p::Vector{<:Real}, basis)
    N  = length(basis)
    dg = zeros(Float64, N, N, N)
    for c in 1:N
        for a in 1:N, b in a:N
            v = dmetric_rotate(p, basis[a], basis[b], basis[c])
            dg[a,b,c] = v
            dg[b,a,c] = v
        end
    end
    return dg
end
