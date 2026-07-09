# src/Metrics.jl
using LinearAlgebra

abstract type AbstractMetric end

"""
    FisherMetric()

Quantum Fisher / Bures metric

    gρ(X,Y) = 1/2 Tr(X Lρ⁻¹(Y))
"""
struct FisherMetric <: AbstractMetric
end

"""
    metric(::FisherMetric, ρ, X, Y)

Evaluate the Fisher/Bures metric

    gρ(X,Y) = 1/2 Tr(X Lρ⁻¹(Y))
"""
function metric(::FisherMetric,
                ρ::AbstractMatrix,
                X::AbstractMatrix,
                Y::AbstractMatrix)

    LY = Lρ_inv(ρ, Y)

    return real(tr(X * LY)) / 2

end

"""
    metric_matrix(metric, ρ, basis)

Construct the Fisher/Bures metric matrix in the given tangent basis.
Optimized via multiple dispatch to separate dense and diagonal density matrices.
"""
function metric_matrix(g::FisherMetric, ρ::AbstractMatrix, basis)
    # Generieke fallback (jouw originele code, maar gecorrigeerd voor symmetrie)
    d = length(basis)
    G = Matrix{Float64}(undef, d, d)
    
    for A in 1:d
        for B in A:d  # Loop alleen over de boven-driehoek
            val = metric(g, ρ, basis[A], basis[B])
            G[A, B] = val
            G[B, A] = val  # Symmetrie uitbuiten
        end
    end
    return G
end

function metric_matrix(g::FisherMetric, ρ::Diagonal{T}, basis) where {T<:Real}
    # Supersnelle route voor diagonaal (bypasst alle SLD matrix constructies)
    d = length(basis)
    mat_dim = size(ρ, 1)
    G = zeros(Float64, d, d)
    diag_ρ = ρ.diag

    # Bereken de inverse denominators één keer vooraf
    inv_denoms = Matrix{T}(undef, mat_dim, mat_dim)
    @inbounds for j in 1:mat_dim
        λ_j = diag_ρ[j]
        for i in 1:mat_dim
            denom = diag_ρ[i] + λ_j
            inv_denoms[i, j] = denom > 1e-12 ? 1.0 / denom : 0.0
        end
    end

    # Element-gewijze contractie zónder SLD allocaties
    for B in 1:d
        TB = basis[B]
        for A in B:d
            TA = basis[A]
            accum = 0.0
            @inbounds for j in 1:mat_dim
                for i in 1:mat_dim
                    accum += real(TA[j, i] * TB[i, j]) * inv_denoms[i, j]
                end
            end
            g_val = 0.5 * accum
            G[A, B] = g_val
            G[B, A] = g_val
        end
    end

    return G
end

function check_metric_normalization(n)

    ρ = maximally_mixed(n)

    basis = gellmann_basis(n)

    g = FisherMetric()

    G = metric_matrix(g, ρ, basis)

    println("‖G - (n/2)I‖ = ",
            norm(G - (n/2)*I, Inf))
end
