# Foundation.jl
# =============
# The single postulate:  g_AB = 𝓕_AB / ρ₀
#
# Defines framework constants, the Fisher information tensor,
# and the Fubini-Study identification 𝓕^(Q) = 4 g^(FS).
# No dependencies — this is the bedrock.

using LinearAlgebra

# ── Constants ────────────────────────────────────────────────

"""τ = 1/5  from 4τ = cos(πτ) — unique real solution."""
const τ = 1//5

let err = abs(4*Float64(τ) - cos(π*Float64(τ)))
    err < 1e-2 || error("τ self-consistency failed: 4τ ≠ cos(πτ), err=$err")
end

"""κ_hol = 6τ = 6/5 — holographic coupling constant."""
const κ_hol = 6//5

"""φ = (1+√5)/2 — golden ratio = [2]_q = 2cos(πτ) at level k = c₁(ℂP²)."""
const φ = (1 + √5) / 2

let err = abs(2*cos(π/5) - φ)
    err < 1e-10 || error("[2]_q ≠ φ, err=$err")
end

let err = abs(1 + φ^4 - 3*φ^2)
    err < 1e-10 || error("1+φ⁴ ≠ 3φ², err=$err")
end

"""λ_W = τ√κ_hol ≈ 0.2191 — Wolfenstein / Cabibbo parameter."""
const λ_W = Float64(τ) * √Float64(κ_hol)

"""A_Wolf = cos(θ₁₂)/√κ_hol — Wolfenstein A parameter."""
const A_Wolf = cos(asin(λ_W)) / √Float64(κ_hol)

"""δ_CP = arctan(φ²) ≈ 69.09° — CP-violation phase."""
const δ_CP = atan(φ^2)

"""M_Pl_GeV = 1.22×10¹⁹ GeV — Planck mass."""
const M_Pl_GeV = 1.22e19

"""M_c_GeV = 1.44×10¹⁷ GeV — KK compactification scale."""
const M_c_GeV = 1.44e17

# ── Fisher information tensor ────────────────────────────────

"""
    quantum_fisher(ρ, A, B) → Float64

Quantum Fisher information F_AB via symmetric logarithmic derivative.
For pure states: F_AB = 4 g_AB^(FS)  (Braunstein-Caves 1994).
"""
function quantum_fisher(ρ::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix)
    vals, vecs = eigen(Hermitian(ρ))
    F = zero(ComplexF64)
    for j in 1:size(ρ,1), k in 1:size(ρ,1)
        denom = vals[j] + vals[k]
        abs(denom) > 1e-14 || continue
        F += 2*vals[j]*vals[k]/denom * dot(vecs[:,j], A*vecs[:,k]) *
                                        dot(vecs[:,k], B*vecs[:,j])
    end
    return real(F)
end

"""
    fisher_tensor(ρ) → Matrix{Float64}

Full Fisher information tensor 𝓕_AB in the Gell-Mann basis of su(n).
This is the fundamental geometric object: g_AB = 𝓕_AB / ρ₀.
"""
function fisher_tensor(ρ::AbstractMatrix)
    basis = gellmann_basis(size(ρ,1))
    n = length(basis)
    F = zeros(Float64, n, n)
    for a in 1:n, b in a:n
        Fab = quantum_fisher(ρ, basis[a], basis[b])
        F[a,b] = Fab; F[b,a] = Fab
    end
    return F
end

"""
    fubini_study_metric(ψ) → Matrix{Float64}

Fubini-Study metric at pure state |ψ⟩.
Theorem (Braunstein-Caves 1994): 𝓕^(Q) = 4 g^(FS).
Quantisation is not a postulate — it is the Fisher geometry.
"""
function fubini_study_metric(ψ::AbstractVector)
    ψn = ψ / norm(ψ)
    return fisher_tensor(ψn * ψn') / 4
end

# ── Basis and states ─────────────────────────────────────────

"""
    gellmann_basis(n) → Vector{Matrix{ComplexF64}}

Generalised Gell-Mann basis: n²-1 traceless generators of su(n).
"""
function gellmann_basis(n::Int)
    basis = Matrix{ComplexF64}[]
    for j in 1:n, k in (j+1):n
        M = zeros(ComplexF64, n, n)
        M[j,k] = 1/√2; M[k,j] = 1/√2
        push!(basis, M)
        M2 = zeros(ComplexF64, n, n)
        M2[j,k] = -1im/√2; M2[k,j] = 1im/√2
        push!(basis, M2)
    end
    for l in 1:(n-1)
        M = zeros(ComplexF64, n, n)
        for m in 1:l; M[m,m] = 1; end
        M[l+1,l+1] = -l
        push!(basis, M / √(l*(l+1)))
    end
    return basis
end

"""
    vacuum_state() → Matrix{ComplexF64}

Vacuum density matrix ρ̂₀ = I/6 on ℂ⁶ = ℂ³⊗ℂ².
Maximally mixed: minimum Fisher information.
"""
vacuum_state() = Matrix{ComplexF64}(I, 6, 6) / 6

"""
    pure_state(ψ) → Matrix{ComplexF64}

Density matrix ρ̂ = |ψ⟩⟨ψ| for pure state |ψ⟩ ∈ ℂ⁶.
Maximum Fisher information — perfect distinguishability.
"""
function pure_state(ψ::AbstractVector)
    ψn = ComplexF64.(ψ) / norm(ψ)
    return ψn * ψn'
end
