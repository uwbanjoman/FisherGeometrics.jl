# lyapunov.jl
# ===========
# Lyapunov-vergelijking voor de Bures-metriek op 𝒟₆
#
# Kern van FisherGeometrics: de metriek g_AB = 𝓕_AB[ρ̂]/ρ₀
# wordt lokaal bepaald door de Bures-metriek via de
# Lyapunov-vergelijking ρ̂G + Gρ̂ = 2X.
#
# Gebruikt door: consciousness.jl, higgs.jl, gravity.jl,
#                en Document CIV (Riemann-hypothese)
#
# © 2026 Jan Bouwman — MIT License

using LinearAlgebra

"""
    lyapunov(ρ, X) → Matrix{ComplexF64}

Oplost de Lyapunov-vergelijking voor de Bures-metriek:

    ρ̂ G + G ρ̂ = 2X  →  G

De Bures-metriek op 𝒟₆ is:

    g_Bures[X, Y] = (1/2) Tr(X G_ρ̂[Y])

waarbij G_ρ̂[Y] = lyapunov(ρ̂, Y).

Algoritme: eigendecompositie van ρ̂, oplossen in diagonaalbasis.

    G_rot[i,j] = 2 X_rot[i,j] / (λ_i + λ_j)

Randgeval: eigenwaarden < 1e-12 worden geregulariseerd.

# Voorbeelden
```julia
ρ = I/6 * ones(6,6)          # vacuüm
G = lyapunov(ρ, G_Y)         # G = 6 G_Y bij vacuüm
g = bures_metric(ρ, G_Y, G_Y) # Bures-norm van G_Y
```

# Verificatie
Bij het vacuüm ρ̂* = I/6 geldt lyapunov(ρ̂*, X) = 6X.
Bewijs: (I/6)G + G(I/6) = G/3 = 2X → G = 6X. ✓
"""
function lyapunov(ρ::AbstractMatrix, X::AbstractMatrix)
    λ, V = eigen(Hermitian(ρ))
    λ = max.(real.(λ), 1e-12)
    X_rot = V' * X * V
    n = size(X, 1)
    G_rot = zeros(ComplexF64, n, n)
    for i in 1:n, j in 1:n
        denom = λ[i] + λ[j]
        G_rot[i,j] = denom > 1e-12 ? 2 * X_rot[i,j] / denom : 0.0
    end
    return V * G_rot * V'
end

"""
    bures_metric(ρ, X, Y) → Float64

Bures-metriek bij ρ̂ in richtingen X, Y:

    g_Bures[X, Y] = (1/2) Tr(X lyapunov(ρ̂, Y))

Dit is de innerproductstructuur op de raakvectorruimte van 𝒟₆
bij het punt ρ̂. Positief definiet voor ρ̂ > 0.

Relatie tot Fisher-informatie:

    g_AB = 𝓕_AB[ρ̂]/ρ₀ = bures_metric(ρ̂, G_A, G_B) / ρ₀
"""
function bures_metric(ρ::AbstractMatrix,
                      X::AbstractMatrix,
                      Y::AbstractMatrix)
    return real(tr(X * lyapunov(ρ, Y))) / 2
end

"""
    bures_norm(ρ, X) → Float64

Bures-norm van richting X bij ρ̂:

    ‖X‖_Bures = √g_Bures[X, X]
"""
function bures_norm(ρ::AbstractMatrix, X::AbstractMatrix)
    return sqrt(max(bures_metric(ρ, X, X), 0.0))
end
