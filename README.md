# FisherGeometrics.jl & M¹·¹·¹-Regularization

A Julia framework for variationally evolving quantum states on the Bures/Fisher
informational manifold toward a geometric vacuum, regularized by the Kaluza-Klein
mass spectrum of the Sasaki-Einstein space M¹·¹·¹.

This repository provides the numerical proof-of-concept that gravity (the macroscopic
Einstein-Hilbert action) can emerge directly from the microscopic quantum information
geometry of states under a variational principle.

---

## 1. Theoretical Architecture

The framework operates on the premise that spacetime geometry is an emergent property
of quantum entanglement and statistical distinguishability.

### 1.1 The Information Basis: The Bures Metric G_ab

Measurable distances and geometry do not arise from physical rulers, but from the
statistical distinguishability of SU(2) density matrices ρ via `FisherMetric()`.
Instead of a predefined spacetime metric, we construct the Bures/Fisher information
metric tensor G directly over SU(2) density matrices ρ using a Gell-Mann tangent basis:

$$G_{ab}(\rho) = 2 \,\text{Re}\!\left[\text{tr}\!\left(\partial_a\rho \cdot \partial_b\rho\right)\right]$$

### 1.2 The Cosmic Regulator: M¹·¹·¹ Sasaki-Einstein Spectroscopy

Quantum fluctuations at the Planck scale introduce wildly divergent infinities.
The framework tames these via the Casimir operator H₀ of the 11-dimensional
Kaluza-Klein compactification on the M¹·¹·¹ space:

$$H_0(M_1, M_2, J) = 4(M_1^2 + M_2^2) + 2J(J+1)$$

This mass term integrates directly into the metric to prevent infinite collapse
(singularities). It yields the **regularized information metric**:

$$G_{\text{reg}} = G_{ab}(\rho) + \alpha H_0 \cdot I$$

High-energy unphysical fluctuations are thereby suppressed by embedding H₀
directly into the metric as a kinetic inertia term.

### 1.3 The Informational Lagrangian

The system dynamics are governed by the informational Einstein-Hilbert Lagrangian
density:

$$\mathcal{L}(\rho, G, \text{Ric}, S) = S(\rho) - 2\Lambda$$

where S(ρ) is the scalar curvature (Ricci scalar) of the information manifold
and Λ is the cosmological constant.

---

## 2. Mathematical Workflow

### 2.1 The Functional Action and Emergent Evolution

The geometry evolves via an action integral analogous to the Einstein-Hilbert action,
driven by the scalar information curvature R(ρ) and the cosmological constant Λ:

$$S = \int_{\mathcal{M}} \bigl(R(\rho) - 2\Lambda\bigr)\,\sqrt{|\det G_{\text{reg}}|}\;d\Omega$$

The actual dynamics — the path toward the stable geometric vacuum — are dictated
by **Quantum Natural Gradient Descent**:

$$\frac{\partial\rho}{\partial t} = -\gamma \cdot G_{\text{reg}}^{-1}(\rho)\,\frac{\delta S}{\delta\rho}$$

The framework drives the state configuration toward the geometric vacuum using
four core mathematical steps:

1. **Numerical Gradient (δS):** Computes the flat variational derivative of
   the total action S = ∫ L dV.

2. **Gell-Mann Projection:** Maps the flat matrix gradient into the covariant
   tangent space of the information manifold via the trace inner product:

   $$h_a = \text{tr}\!\left(\frac{\partial S}{\partial\rho} \cdot \text{basis}_a\right)$$

3. **Natural Gradient & KK-Damping:** Transforms the covariant gradient into
   a contravariant update vector using the inverse of the regularized Fisher metric:

   $$G_{\text{reg}} = G + (\alpha \cdot H_0)\cdot I
   \;\Longrightarrow\;
   h^{\text{contra}} = G_{\text{reg}}^{-1}\cdot h$$

4. **Quantum Projection:** Maps the unconstrained updated matrices back onto the
   valid Bures manifold by enforcing Hermiticity, positivity (eigentruncation),
   and unit trace (tr(ρ) = 1).

### 2.2 At a Glance

| Symbol | Role |
|---|---|
| G_ab(ρ) | Defines the space in which we measure |
| α H₀ | The damper that prevents reality from dissolving into chaos |
| G_reg⁻¹ δS/δρ | Gravity itself: the shortest, most efficient route from a quantum state to its stable geometric vacuum |

---

## 3. Core Package API (`FisherGeometrics`)

The main package exports the following essential functions for the optimization loop:

- **`gellmann_basis(n)`** — Generates the n×n Lie algebra generators.
- **`metric_matrix(g, ρ, basis)`** — Calculates the local 3×3 Fisher metric.
- **`scalar_curvature(g, ρ, basis)`** — Computes the informational Ricci scalar S.
- **`information_action(g, rhos, basis, L; Δ)`** — Integrates the Lagrangian over
  the volume element dV = |det G| Δ.
- **`natural_gradient(g, flat_rhos, basis, δS, rhos_init, M1, M2, J; α)`** —
  Executes the coordinate-aligned, M¹·¹·¹-dampened natural gradient step.
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
    # the 5 iterations may take some seconds to run
    current_action = total_action(flat_current)
    mean_error = sum(abs.(δS)) / length(δS)
    @printf("Iteration %d | Action: %12.6f | Mean Field Error: %10.6f\n", iter, current_action, mean_error)
end

# 4. from the postulate $g_{\mu\nu} = \frac{\mathcal{F}_{\mu\nu}}{\rho_0}$ ​
#    to a Einstein Geometry, induced by the Bures metric
T  = gellmann_basis(6)
FG = bures_einstein(T)

# bures_einstein(T) returns a NamedTuple with:
#    .d      d-symbol tensor (N×N×N)
#    .Γ      Christoffel-symbols (N×N×N)
#    .Q      Riemann ΓΓ-contribution (N×N×N×N)
#    .Ric    Ricci-tensor (N×N)
#    .S_F    Scalar curvature (Float64, expect: 560)
#    .G      Einstein-tensor (N×N)
#    .λ      Ricci-eigenvalue (Float64, expect: 16)
#    .Λ      cosmological constant (Float64, expect: 264)
