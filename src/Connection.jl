# src/Connection.jl

abstract type AbstractConnection end

struct LeviCivitaConnection <: AbstractConnection
end

"""
    christoffel(ρ, basis)

Compute the Levi-Civita connection Γⁱⱼₖ of the Fisher/Bures metric
in the supplied tangent-space basis.

Returns

    Γ[i,j,k] = Γᵏᵢⱼ

where

    ∇_{e_i} e_j = Σₖ Γᵏᵢⱼ e_k
"""
function christoffel(ρ::AbstractMatrix,
                     basis::AbstractVector)

    N = length(basis)

    #
    # metric tensor
    #
    G = metric_matrix(ρ, basis)

    #
    # inverse metric
    #
    Ginv = inv(G)

    #
    # metric derivatives
    #
    D = zeros(Float64,N,N,N)

    g = FisherMetric()

    for k in 1:N
        H = basis[k]

        for i in 1:N
            X = basis[i]

            for j in 1:N
                Y = basis[j]

                D[k,i,j] = dmetric(g,ρ,X,Y,H)

            end
        end
    end

    #
    # Christoffel symbols
    #
    Γ = zeros(Float64,N,N,N)

    for i in 1:N
        for j in 1:N
            for k in 1:N

                s = 0.0

                for l in 1:N

                    s += Ginv[k,l] *
                        (
                            D[i,j,l] +
                            D[j,i,l] -
                            D[l,i,j]
                        )

                end

                Γ[i,j,k] = 0.5*s

            end
        end
    end

    return Γ

end
