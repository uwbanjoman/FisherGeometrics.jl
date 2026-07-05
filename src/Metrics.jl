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

Returns the matrix

    G[A,B] = gρ(T_A, T_B)

where `T_A` and `T_B` are basis elements of the tangent space at `ρ`.
"""
function metric_matrix(::FisherMetric,
                       ρ::AbstractMatrix,
                       basis)

    n = length(basis)

    G = Matrix{Float64}(undef,n,n)

    for A in 1:n
        for B in 1:n
            G[A,B] = metric(FisherMetric(),
                            ρ,
                            basis[A],
                            basis[B])
        end
    end

    return G

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

function check_metric_normalization(n)

    ρ = maximally_mixed(n)

    basis = gellmann_basis(n)

    g = FisherMetric()

    G = metric_matrix(g, ρ, basis)

    println("‖G - (n/2)I‖ = ",
            norm(G - (n/2)*I, Inf))
end
