"""
Evolution.jl
============
Time evolution of the density matrix ρ̂ under the unified equation.

The Von Neumann equation:

    iħ dρ̂/dt = [Ð²_K, ρ̂]

is derived in Document XLIII as the geodesic equation of the
Fisher metric on K — not postulated. The Hamiltonian H = Ð²_K
is uniquely determined by the consistency condition between
the quantum and gravitational sectors (Document XXX).

This module provides:
  - Exact unitary evolution via matrix exponential
  - 4th-order Runge-Kutta integration for general H
  - Information distance and Bures metric along trajectories
  - The measurement as geodesic motion toward maximum 𝓕_AB
  - Entropy and purity tracking
  - Decoherence in the Fisher geometry

The measurement postulate (§9.1 of Document XLIII):
  ρ̂* = argmax 𝓕_AB[ρ̂]

Superposition = low Fisher information (indistinguishable states).
Measurement outcome = high Fisher information (distinguishable state).
Collapse is not an external postulate — it is geodesic motion
toward maximal information distinguishability.

Dependencies: Foundation, Geometry
"""
module Evolution

using LinearAlgebra
using Printf
using ..Foundation: τ, κ_hol, φ, fisher_tensor, quantum_fisher, gellmann_basis
using ..Geometry: bures_distance, spectrum_K

export evolve_exact, evolve_rk4
export von_neumann_rhs
export entropy, purity
export information_distance_trajectory
export decoherence_time
export measurement_projection
export hamiltonian_KK
export check_evolution

# ─────────────────────────────────────────────────────────────
# THE HAMILTONIAN
# ─────────────────────────────────────────────────────────────

"""
    hamiltonian_KK(N) → Matrix{Float64}

Matrix representation of Ð²_K in the truncated KK basis.

The Hamiltonian of the framework is the squared Dirac operator:
    H = Ð²_K

Its eigenvalues are the KK masses M²_n = λ_n + 9/4, where λ_n
are the eigenvalues of the scalar Laplacian Δ_K.

The SM fermions correspond to the lowest eigenspace at λ_min = 9/4.
All Standard Model dynamics follows from the spectrum of this operator.

Arguments:
  N — number of KK levels to include (truncation)

Returns a diagonal N×N matrix with KK mass eigenvalues.
"""
function hamiltonian_KK(N::Int=6)
    λs, _ = spectrum_K(N, N, N)
    M2 = sort(unique(λs .+ 9/4))[1:N]
    return Diagonal(M2)
end

# ─────────────────────────────────────────────────────────────
# TIME EVOLUTION
# ─────────────────────────────────────────────────────────────

"""
    von_neumann_rhs(ρ, H; ħ) → Matrix{ComplexF64}

Right-hand side of the Von Neumann equation:
    dρ̂/dt = −(i/ħ) [H, ρ̂]  =  −(i/ħ)(Hρ̂ − ρ̂H)

This is the geodesic equation of the Fisher metric on K,
not a separate quantum postulate (Document XLIII, Step 3).
"""
function von_neumann_rhs(ρ::AbstractMatrix, H::AbstractMatrix; ħ::Real=1.0)
    return (-1im / ħ) * (H * ρ - ρ * H)
end

"""
    evolve_exact(ρ₀, H, t; ħ) → Matrix{ComplexF64}

Exact unitary time evolution via matrix exponential.

    ρ̂(t) = U(t) ρ̂₀ U†(t)   where   U(t) = exp(−iHt/ħ)

This is the unique solution to the Von Neumann equation for
time-independent H. Unitarity is preserved exactly.

Arguments:
  ρ₀ — initial density matrix
  H  — Hamiltonian (Hermitian matrix)
  t  — evolution time
  ħ  — reduced Planck constant (default 1 in natural units)
"""
function evolve_exact(ρ₀::AbstractMatrix, H::AbstractMatrix, t::Real; ħ::Real=1.0)
    U  = exp(-1im * H * (t / ħ))
    ρt = U * ρ₀ * U'
    # Enforce exact hermiticity (remove floating point drift)
    return (ρt + ρt') / 2
end

"""
    evolve_rk4(ρ₀, H, t_end; dt, ħ) → Vector{Tuple{Float64, Matrix{ComplexF64}}}

4th-order Runge-Kutta integration of the Von Neumann equation.

More flexible than evolve_exact: handles time-dependent H,
open quantum systems (if H is non-Hermitian), and allows
inspection of intermediate states.

Returns a trajectory: vector of (t, ρ̂(t)) pairs.

Arguments:
  ρ₀    — initial density matrix
  H     — Hamiltonian (can be a function H(t) for time-dependent case)
  t_end — final time
  dt    — time step (default 0.01)
  ħ     — reduced Planck constant (default 1)
"""
function evolve_rk4(
    ρ₀::AbstractMatrix,
    H::AbstractMatrix,
    t_end::Real;
    dt::Real=0.01,
    ħ::Real=1.0
)
    trajectory = Tuple{Float64, Matrix{ComplexF64}}[]
    ρ = ComplexF64.(ρ₀)
    t = 0.0
    push!(trajectory, (t, copy(ρ)))

    f(ρ) = von_neumann_rhs(ρ, H; ħ=ħ)

    while t < t_end - dt/2
        k1 = f(ρ)
        k2 = f(ρ + dt/2 * k1)
        k3 = f(ρ + dt/2 * k2)
        k4 = f(ρ + dt   * k3)
        ρ  = ρ + dt/6 * (k1 + 2k2 + 2k3 + k4)
        ρ  = (ρ + ρ') / 2   # enforce hermiticity
        t += dt
        push!(trajectory, (t, copy(ρ)))
    end

    return trajectory
end

# ─────────────────────────────────────────────────────────────
# STATE DIAGNOSTICS
# ─────────────────────────────────────────────────────────────

"""
    purity(ρ) → Float64

Purity of a density matrix: P = Tr(ρ²) ∈ [1/n, 1].

  P = 1    → pure state  (maximum Fisher information)
  P = 1/n  → maximally mixed (minimum Fisher information)

In the framework: purity measures how far ρ̂ is from the
vacuum fixed point ρ̂₀ = I/6 (purity = 1/6 ≈ 0.167).
"""
function purity(ρ::AbstractMatrix)
    return real(tr(ρ * ρ))
end

"""
    entropy(ρ) → Float64

Von Neumann entropy S = −Tr(ρ log ρ).

  S = 0      → pure state
  S = log n  → maximally mixed (vacuum)

The entropy is conserved under the Von Neumann evolution —
unitary evolution cannot create or destroy information.
Information paradoxes arise only when unitarity is violated.
"""
function entropy(ρ::AbstractMatrix)
    vals = real.(eigvals(Hermitian(ρ)))
    return -sum(λ * log(λ) for λ in vals if λ > 1e-15)
end

"""
    information_distance_trajectory(trajectory) → Vector{Float64}

Bures distances along an evolution trajectory.

Returns the distance from the initial state ρ̂₀ to each ρ̂(t):
    d(t) = D_B(ρ̂₀, ρ̂(t))

Useful for visualising how far the state has moved in the
Fisher information geometry during evolution.
"""
function information_distance_trajectory(
    trajectory::Vector{<:Tuple{Float64, Matrix{ComplexF64}}}
)
    ρ₀ = trajectory[1][2]
    return [bures_distance(ρ₀, ρt) for (_, ρt) in trajectory]
end

# ─────────────────────────────────────────────────────────────
# DECOHERENCE
# ─────────────────────────────────────────────────────────────

"""
    decoherence_time(ρ, H, ε; ħ) → Float64

Estimate the decoherence time: how long until off-diagonal
elements of ρ̂ in the H-eigenbasis decay below threshold ε.

For closed quantum systems (unitary evolution), decoherence
does not occur — the off-diagonal elements oscillate but do
not decay. In the FisherGeometrics framework, apparent
decoherence is the approach to a local maximum of 𝓕_AB
on the state space.

Returns Inf for unitary evolution (no true decoherence).
"""
function decoherence_time(
    ρ::AbstractMatrix,
    H::AbstractMatrix,
    ε::Real=0.01;
    ħ::Real=1.0
)
    # For unitary evolution decoherence time is infinite
    # (off-diagonal elements oscillate, not decay)
    vals, vecs = eigen(Hermitian(H))
    ρ_eig = vecs' * ρ * vecs
    max_offdiag = maximum(abs.(ρ_eig - Diagonal(real.(diag(ρ_eig)))))

    if max_offdiag < ε
        return 0.0   # already diagonal in H basis
    else
        return Inf   # unitary evolution: no decay
    end
end

# ─────────────────────────────────────────────────────────────
# MEASUREMENT
# ─────────────────────────────────────────────────────────────

"""
    measurement_projection(ρ, observable) → NamedTuple

Measurement as geodesic motion toward maximum Fisher information.

In the framework (Document XLIII §9.1):
    ρ̂* = argmax 𝓕_AB[ρ̂]

Superposition = low 𝓕_AB = indistinguishable states.
Outcome = high 𝓕_AB = maximally distinguishable.

This function computes:
  1. The eigenvalues of the observable (measurement outcomes)
  2. The projection probabilities P_k = ⟨k|ρ̂|k⟩
  3. The post-measurement states ρ̂_k = |k⟩⟨k|

Arguments:
  ρ          — pre-measurement density matrix
  observable — Hermitian operator being measured

Returns probabilities and post-measurement states.
"""
function measurement_projection(ρ::AbstractMatrix, observable::AbstractMatrix)
    vals, vecs = eigen(Hermitian(observable))
    n = length(vals)

    probabilities = Float64[]
    post_states   = Matrix{ComplexF64}[]

    for k in 1:n
        ψk = vecs[:, k]
        Pk = real(ψk' * ρ * ψk)
        push!(probabilities, max(0.0, Pk))
        push!(post_states, ψk * ψk')
    end

    # Normalise probabilities
    total = sum(probabilities)
    probabilities ./= total

    return (
        outcomes      = vals,
        probabilities = probabilities,
        post_states   = post_states,
    )
end

# ─────────────────────────────────────────────────────────────
# CONSISTENCY CHECKS
# ─────────────────────────────────────────────────────────────

"""
    check_evolution() → Bool

Verify key evolution properties. Returns true if all pass.
"""
function check_evolution()
    ok = true

    # Test 1: Unitarity conservation under exact evolution
    H = [1.0 0.3; 0.3 -1.0]
    ρ₀ = [0.7 0.2+0.1im; 0.2-0.1im 0.3]
    ρ₀ = (ρ₀ + ρ₀') / 2

    ρT = evolve_exact(ρ₀, H, 2π)
    Δtr = abs(tr(ρT) - 1.0)
    if Δtr > 1e-10
        @warn "Tr(ρ) not conserved: deviation $Δtr"
        ok = false
    end

    # Test 2: Purity conservation (unitary evolution)
    P₀ = purity(ρ₀)
    PT = purity(ρT)
    if abs(P₀ - PT) > 1e-10
        @warn "Purity not conserved: $P₀ → $PT"
        ok = false
    end

    # Test 3: Entropy conservation
    S₀ = entropy(ρ₀)
    ST = entropy(ρT)
    if abs(S₀ - ST) > 1e-8
        @warn "Entropy not conserved: $S₀ → $ST"
        ok = false
    end

    # Test 4: Von Neumann RHS is traceless (trace of commutator = 0)
    rhs = von_neumann_rhs(ρ₀, H)
    if abs(tr(rhs)) > 1e-12
        @warn "Tr(dρ/dt) ≠ 0: $(tr(rhs))"
        ok = false
    end

    # Test 5: RK4 agrees with exact evolution
    traj = evolve_rk4(ρ₀, H, 1.0; dt=0.001)
    ρ_rk4 = traj[end][2]
    ρ_exact = evolve_exact(ρ₀, H, 1.0)
    err = maximum(abs.(ρ_rk4 - ρ_exact))
    if err > 1e-6
        @warn "RK4 vs exact deviation: $err"
        ok = false
    end

    return ok
end

end # module Evolution
