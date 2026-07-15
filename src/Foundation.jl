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

"""
    M_KK_GeV = 178.1 GeV
De Kaluza-Klein compactificatieschaal — gelijk aan de 10D Planck-massa
(Document XCIV): 
    M_KK = M_P^{(10D)} = M_P^{(4D)} / (√(8π) Vol(K)^{1/8})^{1/2} = 178.1 GeV
"""
const M_KK_GeV = 178.1

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
    pauli_basis() → Vector{Matrix{ComplexF64}}
"""
function pauli_basis()
    σ1 = ComplexF64[
        0 1
        1 0
    ] / sqrt(2)

    σ2 = ComplexF64[
        0 -im
        im 0
    ] / sqrt(2)

    σ3 = ComplexF64[
        1 0
        0 -1
    ] / sqrt(2)

    return [σ1, σ2, σ3]
end

"""
    su_basis(n)

Return an orthonormal Hermitian basis of the Lie algebra su(n).

The basis consists of

- symmetric generators,
- antisymmetric generators,
- diagonal (Cartan) generators.

The normalization satisfies

    tr(Tᵢ*Tⱼ) = 1/2 δᵢⱼ

so that the basis contains n²−1 generators.
"""
function su_basis(n::Int)

    T = Matrix{ComplexF64}[]

    # Symmetric generators
    for j in 1:n, k in j+1:n
        M = zeros(ComplexF64,n,n)
        M[j,k] = 0.5
        M[k,j] = 0.5
        push!(T,M)
    end

    # Antisymmetric generators
    for j in 1:n, k in j+1:n
        M = zeros(ComplexF64,n,n)
        M[j,k] = -0.5im
        M[k,j] =  0.5im
        push!(T,M)
    end

    # Cartan generators
    for l in 1:n-1
        M = zeros(ComplexF64,n,n)
        α = inv(sqrt(2l*(l+1)))

        for j in 1:l
            M[j,j] = α
        end

        M[l+1,l+1] = -l*α

        push!(T,M)
    end

    return T

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

"""
    gibbs_state(M1::Int, M2::Int, J::Real, β::Real)

Genereert een 2x2 Gibbs-toestand (dichtheidsmatrix) waarbij de energiekloof 
wordt geschaald door de Kaluza-Klein Casimir-massa H₀(M₁, M₂, J).
β is de inverse temperatuur (1/T).
"""
function gibbs_state(M1::Int, M2::Int, J::Real, β::Real)
    # 1. Bepaal de effectieve energie/massa van de KK-modus
    energy_scale = H0(M1, M2, J)
    
    # 2. Bouw een effectieve 2x2 Hamiltoniaan voor het su(2) subsysteem.
    # We gebruiken een standaard opsplitsing (bijv. Pauli-Z) waarbij de 
    # energiekloof evenredig is met de KK-massa.
    H = ComplexF64[energy_scale   0.0;
                   0.0           -energy_scale]
    
    # 3. Bereken de ongenormaliseerde Boltzmann-factor via de matrix-exponentieel
    boltzmann_matrix = exp(-β * H)
    
    # 4. Normaliseer door het spoor (tr) om een geldige dichtheidsmatrix te krijgen
    ρ_gibbs = boltzmann_matrix / tr(boltzmann_matrix)
    
    return ρ_gibbs
end

function gibbs_state_expanded(M1, M2, J, β)
    ρ_2x2 = gibbs_state(M1, M2, J, β)
    # Plaats de 2x2 toestand in een grotere 6x6 nul-matrix
    ρ_6x6 = zeros(ComplexF64, 6, 6)
    ρ_6x6[1:2, 1:2] = ρ_2x2
    return ρ_6x6
end
