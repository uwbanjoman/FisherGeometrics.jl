# src/DensityMatrices.jl

using LinearAlgebra

"""
    maximally_mixed(n)

Return the maximally mixed density matrix I/n.
"""
function maximally_mixed(n::Integer)

    return Matrix{Float64}(I, n, n) / n

end

"""
    is_density_matrix(ρ)

Checks whether ρ is a valid density matrix.
"""
function is_density_matrix(ρ::AbstractMatrix; atol=1e-12)

    # square
    size(ρ,1)==size(ρ,2) || return false

    # Hermitian
    ishermitian(ρ) || return false

    # trace = 1
    abs(tr(ρ)-1) < atol || return false

    # positive semidefinite
    λ = eigvals(Hermitian(ρ))

    minimum(real.(λ)) > -atol || return false

    return true

end
