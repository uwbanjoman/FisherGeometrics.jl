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

function check_metric_normalization(n)

    ρ = maximally_mixed(n)

    basis = gellmann_basis(n)

    g = FisherMetric()

    G = metric_matrix(g, ρ, basis)

    println("‖G - (n/2)I‖ = ",
            norm(G - (n/2)*I, Inf))
end
