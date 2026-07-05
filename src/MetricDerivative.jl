# src/MetricDerivative.jl

using LinearAlgebra

import ..Operators:
    Lρ,
    Lρ_inv,
    dLρ

import ..Metric:
    FisherMetric

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
