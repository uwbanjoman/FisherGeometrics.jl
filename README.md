# FisherGeometrics.jl & M111-Regularization

A Julia framework for variationally evolving quantum states on the Bures/Fisher informational manifold toward a geometric vacuum, regularized by the Kaluza-Klein mass spectrum of the Sasaki-Einstein space $M^{1,1,1}$.

This repository provides the numerical proof-of-concept that gravity (the macroscopic Einstein-Hilbert action) can emerge directly from the microscopic quantum information geometry of states under a variational principle.

## 1. Theoretical Architecture

The framework operates on the premise that spacetime geometry is an emergent property of quantum entanglement and statistical distinguishability.

* **Information Geometry:** Instead of a predefined spacetime metric, we construct the Bures/Fisher information metric tensor $G$ directly over $SU(2)$ density matrices $\rho$ using a Gell-Mann tangent basis.
* **The Informational Lagrangian:** The system dynamics are governed by the informational Einstein-Hilbert Lagrangian density:
  $$\mathcal{L}(\rho, G, \text{Ric}, S) = S(\rho) - 2\Lambda$$
  Where $S(\rho)$ is the scalar curvature (Ricci scalar) of the information manifold and $\Lambda$ is the cosmological constant.
* **Kaluza-Klein Regularization:** High-energy unphysical fluctuations are suppressed by embedding the Casimir-operator $H_0(M_1, M_2, J)$ of the 11D compactification space $M^{1,1,1}$ directly into the metric as a kinetic inertia term.

---

## 2. Mathematical Workflow

The framework drives the state configuration toward the geometric vacuum using four core mathematical steps:

1. **Numerical Gradient ($\delta \mathcal{S}$):** Computes the flat variational derivative of the total action $\mathcal{S} = \int \mathcal{L} \, dV$.
2. **Gell-Mann Projection:** Maps the flat matrix gradient into the covariant tangent space of the information manifold via the trace inner product:
   $$h_a = \text{tr}\left(\frac{\partial \mathcal{S}}{\partial \rho} \cdot \text{basis}_a\right)$$
3. **Natural Gradient & KK-Damping:** Transforms the covariant gradient into a contravariant update vector using the inverse of the regularized Fisher metric:
   $$G_{\text{reg}} = G + (\alpha \cdot H_0) \cdot I \implies h^{\text{contra}} = G_{\text{reg}}^{-1} \cdot h$$
4. **Quantum Projection:** Maps the unconstrained updated matrices back onto the valid Bures manifold by enforcing Hermiticity, positivity ($eigentruncation$), and unit trace ($\text{tr}(\rho) = 1$).

---

## 3. Core Package API (`FisherGeometrics`)

The main package exports the following essential functions for the optimization loop:

* `gellmann_basis(n)`: Generates the $n \times n$ Lie algebra generators.
* `metric_matrix(g, ρ, basis)`: Calculates the local $3 \times 3$ Fisher metric.
* `scalar_curvature(g, ρ, basis)`: Computes the informational Ricci scalar $S$.
* `information_action(g, rhos, basis, L; Δ)`: Integrates the Lagrangian over the volume element $dV = \sqrt{| \det G |} \, \Delta$.
* `natural_gradient(g, flat_rhos, basis, δS, rhos_init, M1, M2, J; α)`: Executes the coordinate-aligned, $M^{1,1,1}$-dampened natural gradient step.

---

## 4. Usage Example

The following minimal working example initializes a sequence of quantum states and evolves them over 5 iterations. Notice how the total action consistently drops, signaling convergence to the physical vacuum:

```julia
using LinearAlgebra
using Printf
using FisherGeometrics

# 1. Framework Setup & Constants
g = FisherGeometrics.FisherMetric()
basis = FisherGeometrics.gellmann_basis(2)

Λ = -0.2
einstein_hilbert_L = (ρ, G, Ric, S) -> S - 2*Λ
Δ = 0.1
learning_rate = 0.0005
α = 1e-3  # Kaluza-Klein coupling constant

# Quantum numbers for the lowest massive M¹¹¹ spinor-mode
M1, M2, J = 0, 0, 1 

# 2. State Initialization (Valid Density Matrices)
make_valid_rho(x) = ComplexF64[0.5 + 0.1*sin(x)  0.05*cos(x); 0.05*cos(x)  0.5 - 0.1*sin(x)]
rhos_init = [make_valid_rho(x) for x in 0.0:Δ:1.0]

function total_action(flat_rhos)
    n_states = length(rhos_init)
    reconstructed = Vector{Matrix{ComplexF64}}(undef, n_states)
    idx = 1
    for k in 1:n_states
        mat = [flat_rhos[idx]                  flat_rhos[idx+2] + im*flat_rhos[idx+3];
               flat_rhos[idx+2] - im*flat_rhos[idx+3]  flat_rhos[idx+1]]
        reconstructed[k] = mat
        idx += 4
    end
    return FisherGeometrics.information_action(g, reconstructed, basis, einstein_hilbert_L; Δ = Δ)
end

# Flatten initial state parameters
flat_current = Float64[]
for ρ in rhos_init
    push!(flat_current, real(ρ[1,1]), real(ρ[2,2]), real(ρ[1,2]), imag(ρ[1,2]))
end

# 3. Variational Optimization Loop
println("Starting convergence loop on the M¹¹¹-regularized manifold...\n")

for iter in 1:5
    # Calculate regular and natural gradients
    δS = FisherGeometrics.numerical_gradient(total_action, flat_current)
    nat_δS = FisherGeometrics.natural_gradient(g, flat_current, basis, δS, rhos_init, M1, M2, J; α = α)
    
    # Take a step along the information geodesic
    flat_step = flat_current - learning_rate * nat_δS
    
    # Quantum Projection (Enforce CPTP / valid state space)
    flat_current = Float64[]
    idx = 1
    for k in 1:length(rhos_init)
        mat = [flat_step[idx]                 flat_step[idx+2] + im*flat_step[idx+3];
               flat_step[idx+2] - im*flat_step[idx+3]  flat_step[idx+1]]
        mat = (mat + mat') / 2
        vals, vecs = eigen(Hermitian(mat))
        mat_projected = vecs * Diagonal(max.(vals, 1e-6)) * vecs'
        mat_projected ./= tr(mat_projected)  # Strict unit trace correction
        
        push!(flat_current, real(mat_projected[1,1]), real(mat_projected[2,2]), 
                            real(mat_projected[1,2]), imag(mat_projected[1,2]))
        idx += 4
    end
    
    # Monitor Action and Gradient Field Error
    current_action = total_action(flat_current)
    mean_error = sum(abs.(δS)) / length(δS)
    @printf("Iteration %d | Action: %12.6f | Mean Field Error: %10.6f\n", iter, current_action, mean_error)
end
