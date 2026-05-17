# Evolution.jl
# ============
# Von Neumann time evolution: iħ dρ̂/dt = [Ð²_K, ρ̂]
# Derived as geodesic equation of the Fisher metric (Document XLIII).
# Depends on: Foundation.jl, Geometry.jl

# ── Hamiltonian ───────────────────────────────────────────────

"""
    hamiltonian_KK(N) → Diagonal{Float64}

Ð²_K in truncated KK basis: diagonal with eigenvalues M²_n = λ_n + 9/4.
SM fermions at lowest eigenvalue 9/4.
"""
function hamiltonian_KK(N::Int=6)
    M2 = kk_masses(N)
    return Diagonal(M2)
end

# ── Time evolution ────────────────────────────────────────────

"""
    von_neumann_rhs(ρ, H; ħ) → Matrix{ComplexF64}

dρ̂/dt = −(i/ħ)[H, ρ̂]. Geodesic equation of Fisher metric — not a postulate.
"""
von_neumann_rhs(ρ::AbstractMatrix, H::AbstractMatrix; ħ::Real=1.0) =
    (-1im/ħ) * (H*ρ - ρ*H)

"""
    evolve_exact(ρ₀, H, t; ħ) → Matrix{ComplexF64}

Exact unitary evolution: ρ̂(t) = e^{-iHt/ħ} ρ̂₀ e^{+iHt/ħ}.
"""
function evolve_exact(ρ₀::AbstractMatrix, H::AbstractMatrix, t::Real; ħ::Real=1.0)
    U  = exp(-1im * H * (t/ħ))
    ρt = U * ρ₀ * U'
    return (ρt + ρt') / 2
end

"""
    evolve_rk4(ρ₀, H, t_end; dt, ħ) → Vector{Tuple{Float64, Matrix{ComplexF64}}}

4th-order Runge-Kutta integration. Returns trajectory [(t, ρ̂(t)), ...].
"""
function evolve_rk4(ρ₀::AbstractMatrix, H::AbstractMatrix, t_end::Real;
                    dt::Real=0.01, ħ::Real=1.0)
    traj = Tuple{Float64,Matrix{ComplexF64}}[]
    ρ = ComplexF64.(ρ₀); t = 0.0
    push!(traj, (t, copy(ρ)))
    f(ρ) = von_neumann_rhs(ρ, H; ħ)
    while t < t_end - dt/2
        k1=f(ρ); k2=f(ρ+dt/2*k1); k3=f(ρ+dt/2*k2); k4=f(ρ+dt*k3)
        ρ = (ρ + dt/6*(k1+2k2+2k3+k4)); ρ = (ρ+ρ')/2
        t += dt; push!(traj, (t, copy(ρ)))
    end
    return traj
end

# ── State diagnostics ─────────────────────────────────────────

"""
    purity(ρ) → Float64

Tr(ρ²) ∈ [1/n, 1]. Purity=1: pure state (max Fisher info). Purity=1/n: vacuum.
"""
purity(ρ::AbstractMatrix) = real(tr(ρ*ρ))

"""
    entropy(ρ) → Float64

Von Neumann entropy S = −Tr(ρ log ρ). Conserved under unitary evolution.
"""
function entropy(ρ::AbstractMatrix)
    vals = real.(eigvals(Hermitian(ρ)))
    return -sum(λ*log(λ) for λ in vals if λ > 1e-15)
end

"""
    information_distance_trajectory(traj) → Vector{Float64}

Bures distances from ρ̂₀ to each ρ̂(t) along trajectory.
"""
function information_distance_trajectory(traj::Vector{<:Tuple{Float64,Matrix{ComplexF64}}})
    ρ₀ = traj[1][2]
    return [bures_distance(ρ₀, ρt) for (_, ρt) in traj]
end

# ── Measurement ───────────────────────────────────────────────

"""
    measurement_projection(ρ, observable) → NamedTuple

Measurement as projection onto eigenstates of observable.
Probabilities P_k = ⟨k|ρ̂|k⟩, post-states |k⟩⟨k|.
In the framework: collapse = geodesic motion toward max 𝓕_AB.
"""
function measurement_projection(ρ::AbstractMatrix, observable::AbstractMatrix)
    vals, vecs = eigen(Hermitian(observable))
    probs = [max(0.0, real(vecs[:,k]'*ρ*vecs[:,k])) for k in 1:length(vals)]
    probs ./= sum(probs)
    return (outcomes=vals, probabilities=probs,
            post_states=[vecs[:,k]*vecs[:,k]' for k in 1:length(vals)])
end

"""
    decoherence_time(ρ, H, ε; ħ) → Float64

For unitary evolution returns Inf (no true decoherence in closed system).
"""
function decoherence_time(ρ::AbstractMatrix, H::AbstractMatrix,
                          ε::Real=0.01; ħ::Real=1.0)
    vals, vecs = eigen(Hermitian(H))
    ρ_eig = vecs'*ρ*vecs
    max_off = maximum(abs.(ρ_eig - Diagonal(real.(diag(ρ_eig)))))
    return max_off < ε ? 0.0 : Inf
end

# ── Consistency check ─────────────────────────────────────────

function check_evolution()
    ok = true
    H  = [1.0 0.3; 0.3 -1.0]
    ρ₀ = let r=[0.7 0.2+0.1im; 0.2-0.1im 0.3]; (r+r')/2; end

    abs(tr(von_neumann_rhs(ρ₀,H))) < 1e-12 ||
        (@warn "Tr(dρ/dt) ≠ 0"; ok = false)

    ρT = evolve_exact(ρ₀, H, 2π)
    abs(tr(ρT)-1.0) < 1e-10   || (@warn "Tr(ρ) not conserved"; ok = false)
    abs(purity(ρ₀)-purity(ρT)) < 1e-10 || (@warn "Purity not conserved"; ok = false)
    abs(entropy(ρ₀)-entropy(ρT)) < 1e-8 || (@warn "Entropy not conserved"; ok = false)

    traj = evolve_rk4(ρ₀, H, 1.0; dt=0.001)
    ρ_rk4  = traj[end][2]
    ρ_exact = evolve_exact(ρ₀, H, 1.0)
    maximum(abs.(ρ_rk4-ρ_exact)) < 1e-6 ||
        (@warn "RK4 vs exact deviation > 1e-6"; ok = false)

    return ok
end
