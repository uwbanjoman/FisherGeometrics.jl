---
title: Getting Started
nav_order: 2
---

# Getting Started

## Installation

```julia
] add https://github.com/uwbanjoman/FisherGeometrics.jl
```

## Run the tests

```julia
] test FisherGeometrics
# 89 tests pass
```

## First steps

```julia
using FisherGeometrics

# The vacuum state ρ̂* = I/6
ρ_vac = vacuum_state()

# A pure state
ψ = ComplexF64[1, 0, 0, 0, 0, 0]
ρ = pure_state(ψ)

# Fisher information
fisher_tensor(ρ)           # 35×35 tensor 𝓕_AB
fisher_excess(ρ)           # 5/6 for pure states

# Bures distance
bures_distance(ρ, ρ_vac)  # arccos(1/√6) ≈ 65.9° for any pure state

# Von Neumann evolution
H = Matrix(hamiltonian_KK(6))
ρ_t = evolve_exact(ρ, H, 1.0)
```

## Run an example

```bash
julia --project examples/yang_mills.jl
julia --project examples/black_holes.jl
julia --project examples/gravitational_waves.jl
```

## The single postulate

Everything in this framework follows from:

$$g_{AB} = \frac{\mathcal{F}_{AB}}{\rho_0}$$

The spacetime metric **is** the quantum Fisher information tensor.
See [The Postulate](framework/postulate) for details.
