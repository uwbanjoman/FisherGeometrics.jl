# src/DensityMatrices.jl

using LinearAlgebra

"""
    maximally_mixed(n)

Return the maximally mixed density matrix I/n.
"""
function maximally_mixed(n::Integer)

    return Matrix{Float64}(I, n, n) / n

end
