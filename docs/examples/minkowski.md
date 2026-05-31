---
title: Minkowski Spacetime
parent: Examples
nav_order: 1
---

# Minkowski Spacetime from Conscious Agents

Derives η_μν = diag(−1,+1,+1,+1) from the Fisher information
geometry of N conscious agents in the classical limit.

## Run it

```bash
julia --project examples/spacetime_simulation.jl
```

## What it demonstrates

1. **Covariance metric** — Cov(G_μ, G_ν) = (1/3)δ_μν in vacuum
2. **Wick rotation** — G_t is the U(1)_Y singlet (Y=0) → time gets −
3. **Minkowski** — η_μν = diag(−1,+1,+1,+1)

## The key result

All pure states are equidistant from the vacuum:

$$D_B(|\psi\rangle\langle\psi|,\, \mathbf{I}/6) = \arccos(1/\sqrt{6})
\approx 65.9°$$

for any normalised ψ ∈ ℂ⁶. This is the geometric origin of the
universality of physics — all observers start equally far from vacuum.

## Code

```julia
using FisherGeometrics

# Vacuum state
ρ = vacuum_state()

# Covariance metric → (1/3)δ_μν
G = make_generators()
g = [real(tr(ρ*G[i]*G[j])) - real(tr(ρ*G[i]))*real(tr(ρ*G[j]))
     for i in 1:4, j in 1:4]

# After Wick rotation and normalisation: η_μν = diag(-1,+1,+1,+1)
```

## Related documents

- Document LXXIX — Spacetime from Conscious Agents
- Document LXXX — Minkowski from N Agents (quantum CLT)
