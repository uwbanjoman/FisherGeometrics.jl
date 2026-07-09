# src/Connection.jl

abstract type AbstractConnection end

struct LeviCivitaConnection <: AbstractConnection
end

"""
    christoffel(g, ρ, basis)

Compute the Levi-Civita connection Γⁱⱼₖ of the Fisher/Bures metric
in the supplied tangent-space basis.

Returns

    Γ[i,j,k] = Γᵏᵢⱼ

where

    ∇_{e_i} e_j = Σₖ Γᵏᵢⱼ e_k
"""
function christoffel(g::FisherMetric, ρ, basis)

    G    = metric_matrix(g, ρ, basis)
    Ginv = pinv(G)
    dg   = metric_derivatives(g, ρ, basis)

    nbasis = length(basis)

    Γ = zeros(eltype(G), nbasis, nbasis, nbasis)

    for i in 1:nbasis
        for j in 1:nbasis
            for k in 1:nbasis
                Γ[i,j,k] = 0.5 * sum(
                    Ginv[k,l] *
                    (dg[j,l,i] + dg[i,l,j] - dg[i,j,l])
                    for l in 1:nbasis
                )
            end
        end
    end

    return Γ
end
