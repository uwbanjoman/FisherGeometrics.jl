# src/Curvature.jl

"""
    riemann(g, ρ, basis)

Compute the Riemann curvature tensor

    R[m,i,j,k] = Rᵐᵢⱼₖ

associated with the Fisher metric.

Returns an n×n×n×n array, where
n = length(basis).
"""
function riemann(g::FisherMetric, ρ::AbstractMatrix, basis)
    n = length(basis)

    G    = metric_matrix(g, ρ, basis)
    Ginv = pinv(G)

    dg   = metric_derivatives(g, ρ, basis)   # n × n × n
    ddg  = ddmetric_tensor(g, ρ, basis)     # n × n × n × n
    Γ    = christoffel(g, ρ, basis)          # n × n × n

    # --- OPTIMALISATIE 1: Bereken dGinv vooraf via matrix-multiplicaties ---
    # dGinv_tensor[m, l, i] = ∂ᵢ G⁻¹_{ml}
    dGinv_tensor = zeros(Float64, n, n, n)
    for i in 1:n
        # Pak de slice van dg in de i-richting: dg[:, :, i]
        dg_i = @view dg[:, :, i]
        # De contractie -Ginv * dg_i * Ginv via geoptimaliseerde matrix-multiplicatie
        dGinv_tensor[:, :, i] .= -Ginv * dg_i * Ginv
    end

    dΓ = zeros(Float64, n, n, n, n)

    # --- OPTIMALISATIE 2: Gereduceerde 5D-loop ipv 7D-loop ---
    for k in 1:n
        for j in 1:n
            # Handige shortcuts om herhaaldelijke array-lookups in de loops te vermijden
            dg_k_j = @view dg[k, :, j]
            dg_j_k = @view dg[j, :, k]
            dg_j_k_l = @view dg[j, k, :]
            
            for i in 1:n
                ddg_k_j_i = @view ddg[k, :, j, i]
                ddg_j_k_i = @view ddg[j, :, k, i]
                ddg_j_k_l_i = @view ddg[j, k, :, i]
                
                for m in 1:n
                    s = 0.0
                    @simd for l in 1:n
                        # Gebruik de vooraf berekende dGinv_tensor (O(1) lookup!)
                        term1 = dGinv_tensor[m, l, i] * (dg_k_j[l] + dg_j_k[l] - dg_j_k_l[l])
                        
                        term2 = Ginv[m, l] * (ddg_k_j_i[l] + ddg_j_k_i[l] - ddg_j_k_l_i[l])
                        
                        s += term1 + term2
                    end
                    dΓ[m, i, j, k] = 0.5 * s
                end
            end
        end
    end

    # --- Riemann Tensor Assemblage ---
    R = zeros(Float64, n, n, n, n)
    for k in 1:n
        for j in 1:n
            for i in 1:n
                for m in 1:n
                    val = dΓ[m, i, j, k] - dΓ[m, j, i, k]
                    
                    # Contractie van de Christoffel-symbolen
                    s_Γ = 0.0
                    @simd for l in 1:n
                        s_Γ += Γ[m, i, l] * Γ[l, j, k] - Γ[m, j, l] * Γ[l, i, k]
                    end
                    
                    R[m, i, j, k] = val + s_Γ
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
function scalar_curvature(g::FisherMetric, ρ::AbstractMatrix, basis)
    if ρ isa Diagonal
        # Als het al een Diagonal is, mag hij HIER niet inkomen.
        # Dit vangt situaties op waarin de splitsing misgaat.
    elseif isdiag(ρ)
        return scalar_curvature(g, Diagonal(ρ), basis)
    end
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
