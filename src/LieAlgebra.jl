# LieAlgebra.jl

function commutator(A::Matrix{ComplexF64},B::Matrix{ComplexF64})
  A*B - B*A
end

function anticommutator(A::Matrix{ComplexF64},B::Matrix{ComplexF64})
  A*B + B*A
end

function inner(A::Matrix{ComplexF64},B::Matrix{ComplexF64})
  real(tr(A*B))
end

function structure_constants(G::Vector{Matrix{ComplexF64}})
    n = length(G)
    f = zeros(Float64, n, n, n)
    for a in 1:n
        for b in 1:n
            # commutator
            C = G[a]*G[b] - G[b]*G[a]
            for c in 1:n
                # projecteer op generator c
                #f[a,b,c] = -imag(tr(C * G[c])) / sqrt(2)
                f[a,b,c] = imag(tr(C * G[c]))  # + teken, geen 1/√2
                #f[a,b,c] = -imag(tr(C * G[c]))
            end
        end
    end
    return f
end

function jacobi_test(f::Array{Float64,3}; atol=1e-12)
    n = size(f,1)
    max_error = 0.0
    for a in 1:n
        for b in 1:n
            for c in 1:n
                for e in 1:n
                    s = 0.0
                    for d in 1:n
                        s +=
                            f[a,b,d]*f[d,c,e] +
                            f[b,c,d]*f[d,a,e] +
                            f[c,a,d]*f[d,b,e]
                    end
                    max_error = max(max_error, abs(s))
                end
            end
        end
    end

    println("Maximum Jacobi error = ", max_error)

    return max_error < atol

end

function reconstruction_error(G, f)
    n = length(G)
    maxerr = 0.0
    for a in 1:n
        for b in 1:n
            lhs = commutator(G[a], G[b])
            rhs = zeros(ComplexF64, size(lhs))
            for c in 1:n
                #rhs += im * sqrt(2) * f[a,b,c] * G[c]
                rhs += im * f[a,b,c] * G[c]
                #rhs += im * f[a,b,c] * G[c]
            end
            maxerr = max(maxerr, norm(lhs - rhs))
        end
    end
    return maxerr
end

function expand_in_basis(X, G)
    return [real(tr(X*Gi)) for Gi in G]
end
