# src/Geodesics.jl

"""
    geodesic_equation(g, ρ, basis, v)

Returns the acceleration

    aᵏ = -Γᵏᵢⱼ vⁱ vʲ

where `v` are the coordinates of the tangent vector in the chosen basis.
"""
function geodesic_equation(g::FisherMetric,
                           ρ::AbstractMatrix,
                           basis,
                           v::AbstractVector)

    Γ = christoffel(g, ρ, basis)

    n = length(basis)

    a = zeros(Float64,n)

    for k in 1:n

        s = 0.0

        for i in 1:n
            for j in 1:n

                s += Γ[k,i,j] * v[i] * v[j]

            end
        end

        a[k] = -s

    end

    return a

end
