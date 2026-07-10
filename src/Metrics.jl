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

"""
    natural_gradient(g, flat_rhos, basis, δS, rhos_init, M1, M2, J; α=1e-4)

Berekent de Natural Gradient op de Bures-variëteit, geregulariseerd door 
de Kaluza-Klein Casimir-operator H₀ van M^{1,1,1}.

α is de koppelingsconstante die bepaalt hoe sterk de KK-massa de rimpelingen onderdrukt.
"""
function natural_gradient(g::FisherGeometrics.FisherMetric, flat_rhos, basis, δS, rhos_init, M1::Int, M2::Int, J::Real; α=1e-4)
    n_states = length(rhos_init)
    dim = 2
    h0_weight = H0(M1, M2, J)
    
    # 1. Herbouw de actuele toestanden
    reconstructed = Vector{Matrix{ComplexF64}}(undef, n_states)
    idx = 1
    for k in 1:n_states
        mat = zeros(ComplexF64, dim, dim)
        mat[1,1] = flat_rhos[idx]
        mat[2,2] = flat_rhos[idx+1]
        mat[1,2] = flat_rhos[idx+2] + im*flat_rhos[idx+3]
        mat[2,1] = flat_rhos[idx+2] - im*flat_rhos[idx+3]
        reconstructed[k] = mat
        idx += 4
    end
    
    nat_grad = zero(δS)
    idx = 1
    for k in 1:n_states
        ρ = reconstructed[k]
        
        # Bereken de echte 3x3 Fisher-metriek in de Gell-Mann-basis
        G = metric_matrix(g, ρ, basis)
        G_reg = G + (α * h0_weight) * I
        Ginv = pinv(G_reg)
        
        # Haal de 4 matrix-gradiënt componenten op
        δS_k = δS[idx:(idx+3)]
        
        # Bouw de gradiënt-matrix dS/dρ
        dS_dρ = [δS_k[1]                 δS_k[3] + im*δS_k[4];
                 δS_k[3] - im*δS_k[4]    δS_k[2]]
        
        # PROJECTIE 1: Projecteer dS/dρ naar de covariante Gell-Mann componenten (hulpvector h)
        # h_a = tr(dS_dρ' * basis_a)
        h = zeros(3)
        for a in 1:3
            h[a] = real(tr(dS_dρ * basis[a])) 
        end
        
        # METRIEKE CORRECTIE: Haal de index op in de Gell-Mann-ruimte
        h_contravariant = Ginv * h
        
        # PROJECTIE 2: Transformeer de gecorrigeerde vector terug naar een matrix-richting
        dρ_nat = zeros(ComplexF64, dim, dim)
        for a in 1:3
            dρ_nat += h_contravariant[a] * basis[a]
        end
        
        # Sla de gecorrigeerde componenten op in de vlakke output-vector
        nat_grad_k = zeros(4)
        nat_grad_k[1] = real(dρ_nat[1,1])
        nat_grad_k[2] = real(dρ_nat[2,2])
        nat_grad_k[3] = real(dρ_nat[1,2])
        nat_grad_k[4] = imag(dρ_nat[1,2])
        
        nat_grad[idx:(idx+3)] .= nat_grad_k
        idx += 4
    end
    return nat_grad
end

function check_metric_normalization(n)

    ρ = maximally_mixed(n)

    basis = gellmann_basis(n)

    g = FisherMetric()

    G = metric_matrix(g, ρ, basis)

    println("‖G - (n/2)I‖ = ",
            norm(G - (n/2)*I, Inf))
end
