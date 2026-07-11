# src/Christoffel.jl

function christoffel_symbols(g::FisherMetric, ρ::Matrix, basis::Vector)
    n = length(basis)
    # 1. Bereken de metriek op punt ρ
    G = metric_matrix(g, ρ, basis)
    G_inv = inv(G)
    
    # 2. Haal de partiële afgeleiden van de metriek op (n x n x n tensor)
    # ∂g[m, i, j] = ∂g_ij / ∂ρ_m
    ∂g = metric_derivatives(g, ρ, basis) 
    
    # 3. Initialiseer de Christoffel tensor
    Γ = zeros(n, n, n)
    
    # 4. Berekening volgens de metriek-compatibele connectie
    for k in 1:n
        for i in 1:n
            for j in 1:n
                sum_term = 0.0
                for m in 1:n
                    # De formule: 1/2 * G_inv * (∂_i g_mj + ∂_j g_mi - ∂_m g_ij)
                    sum_term += G_inv[m, k] * (∂g[i, m, j] + ∂g[j, m, i] - ∂g[m, i, j])
                end
                Γ[k, i, j] = 0.5 * sum_term
            end
        end
    end
    
    return Γ
end
