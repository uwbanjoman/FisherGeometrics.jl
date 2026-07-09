# src/Curvature.jl

"""
    riemann(g, ρ, basis)

Compute the Riemann curvature tensor

    R[m,i,j,k] = Rᵐᵢⱼₖ

associated with the Fisher metric.

Returns an n×n×n×n array, where
n = length(basis).
"""
function riemann(g::FisherMetric,
                 ρ::AbstractMatrix,
                 basis)

    n = length(basis)

    G    = metric_matrix(g, ρ, basis)
    Ginv = pinv(G)

    dg  = metric_derivatives(g, ρ, basis)
    ddg = ddmetric_tensor(g, ρ, basis)

    Γ = christoffel(g, ρ, basis)

    dΓ = zeros(Float64,n,n,n,n)

    #
    # dΓ[m,i,j,k]
    #
    # = ∂ᵢ Γᵐⱼₖ
    #

    for m in 1:n
        for i in 1:n
            for j in 1:n
                for k in 1:n

                    s = 0.0

                    for l in 1:n

                        #
                        # ∂ᵢ G^{-1}
                        #

                        dGinv =
                            0.0

                        for a in 1:n
                            for b in 1:n

                                dGinv -=
                                    Ginv[m,a] *
                                    dg[a,b,i] *
                                    Ginv[b,l]

                            end
                        end

                        term1 =
                            dGinv *
                            (
                                dg[k,l,j] +
                                dg[j,l,k] -
                                dg[j,k,l]
                            )

                        term2 =
                            Ginv[m,l] *
                            (
                                ddg[k,l,j,i] +
                                ddg[j,l,k,i] -
                                ddg[j,k,l,i]
                            )

                        s += 0.5*(term1 + term2)

                    end

                    dΓ[m,i,j,k] = s

                end
            end
        end
    end

    #
    # Riemann tensor
    #

    R = zeros(Float64,n,n,n,n)

    for m in 1:n
        for i in 1:n
            for j in 1:n
                for k in 1:n

                    R[m,i,j,k] =
                        dΓ[m,i,j,k] -
                        dΓ[m,j,i,k]

                    for l in 1:n

                        R[m,i,j,k] +=
                            Γ[m,i,l]*Γ[l,j,k] -
                            Γ[m,j,l]*Γ[l,i,k]

                    end

                end
            end
        end
    end

    return R

end

"""
    ricci(g, ρ, basis)

Compute the Ricci curvature tensor associated with the
Fisher metric.

Returns an n×n matrix

    Ric[i,j].
"""
function ricci(g::FisherMetric,
               ρ::AbstractMatrix,
               basis)

    R = riemann(g, ρ, basis)

    n = length(basis)

    Ric = zeros(Float64,n,n)

    for i in 1:n
        for j in 1:n

            s = 0.0

            for m in 1:n
                s += R[m,i,m,j]
            end

            Ric[i,j] = s

        end
    end

    return Ric

end

"""
    scalar_curvature(g, ρ, basis)

Fisher scalar curvature

S = g^{ij} R_{ij}
"""
function scalar_curvature(g::FisherMetric,
                          ρ::AbstractMatrix,
                          basis)

    G    = metric_matrix(g, ρ, basis)
    Ginv = pinv(G)

    Ric = ricci(g, ρ, basis)

    n = size(G,1)

    S = 0.0

    for i in 1:n
        for j in 1:n
            S += Ginv[i,j] * Ric[i,j]
        end
    end

    return S

end
