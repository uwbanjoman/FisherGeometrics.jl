# src/MetricDerivative.jl

using LinearAlgebra

"""
    dmetric(metric, ρ, X, Y, Z)

Directional derivative of the Fisher/Bures metric
along the tangent direction `Z`.

Computes

    D_Z gρ(X,Y)

using

    D(Lρ⁻¹)
        =
        -Lρ⁻¹ ∘ D(Lρ) ∘ Lρ⁻¹.
"""
function dmetric(::FisherMetric,
                 ρ::AbstractMatrix,
                 X::AbstractMatrix,
                 Y::AbstractMatrix,
                 Z::AbstractMatrix)

    # Lρ⁻¹(Y)
    LY = Lρ_inv(ρ, Y)

    # D(Lρ)(LY)
    DLY = dLρ(Z, LY)

    # Lρ⁻¹(D(Lρ)(LY))
    term = Lρ_inv(ρ, DLY)

    # minus sign from
    #
    # D(L⁻¹)=−L⁻¹(DL)L⁻¹
    #
    return -real(tr(X * term))/2

end

"""
    metric_derivatives(g, ρ, basis)

Compute the directional derivatives of the Fisher metric in the
basis `basis`.

Returns a 3-tensor `dg` with

    dg[i,j,k] = Dgρ(H_k; X_i, X_j),

where

- `X_i = basis[i]`,
- `H_k = basis[k]`.

The tensor `dg` is symmetric in its first two indices.
"""
function metric_derivatives(g::FisherMetric,
                            ρ::AbstractMatrix,
                            basis)

    n = length(basis)

    dg = zeros(Float64, n, n, n)

    for i in 1:n
        Xi = basis[i]

        for j in i:n
            Xj = basis[j]

            for k in 1:n
                Hk = basis[k]

                val = dmetric(g, ρ, Xi, Xj, Hk)

                dg[i,j,k] = val
                dg[j,i,k] = val
            end
        end
    end

    return dg
end

"""
    ddmetric(g, ρ, X, Y, H, K)

Second directional derivative of the Fisher metric.

Returns

    D²gρ(H,K)(X,Y)

where `H` and `K` are tangent directions and `X`,`Y` are the metric
arguments.

The implementation uses

    D²Lρ = 0

since the Jordan operator is linear in `ρ`.
"""
function ddmetric(::FisherMetric,
                  ρ::AbstractMatrix,
                  X::AbstractMatrix,
                  Y::AbstractMatrix,
                  H::AbstractMatrix,
                  K::AbstractMatrix)

    # Lρ⁻¹(Y)
    LY = Lρ_inv(ρ, Y)

    #
    # First contribution
    #
    T1 = Lρ_inv(ρ,
                jordan(K,
                    Lρ_inv(ρ,
                        jordan(H, LY))))

    #
    # Second contribution
    #
    T2 = Lρ_inv(ρ,
                jordan(H,
                    Lρ_inv(ρ,
                        jordan(K, LY))))

    return real(tr(X * (T1 + T2))) / 2

end
