---
title: Home
nav_order: 1
---

# FisherGeometrics.jl

**One density matrix. One postulate. Everything.**

$$g_{AB} = \mathcal{F}_{AB} / \rho_0$$

The spacetime metric **is** the quantum Fisher information tensor.
All Standard Model structure follows as a geometric consequence.
Zero free parameters.

---

## What this framework derives

| Result | Value | Deviation |
|--------|-------|-----------|
| Weinberg angle sin²θ_W | 0.232 | 0.3% |
| Fine structure constant 1/α | 137.08 | 0.03% |
| CP phase δ_CP | 69.09° | 0.15% |
| Three generations | 3 | exact |
| Hypercharges | +1/6, +2/3, −1/3 | exact |
| Yang-Mills mass gap Δ | 9/4 | exact |
| Bekenstein-Hawking entropy | A/4G_N | exact |
| Neutrino mass sum Σmν | 58.7 meV | testable |

All from a single postulate. Zero free parameters.

---

## Falsifiable prediction

$$\Sigma m_\nu \approx 58.7 \text{ meV}$$

Testable by Euclid (2025–2030).
A measurement of Σmν < 30 meV falsifies the framework.

---

## Quick start

```julia
] add https://github.com/uwbanjoman/FisherGeometrics.jl

using FisherGeometrics

# The fundamental object
ρ = pure_state([1,0,0,0,0,0])

# Fisher information
F = fisher_tensor(ρ)

# Bures distance
d = bures_distance(ρ, vacuum_state())

# Von Neumann evolution
H = hamiltonian_KK(6)
ρ_t = evolve_exact(ρ, H, 1.0)
```

---

## Navigation

- [Getting started](getting_started) — installation and first steps
- [The framework](framework/postulate) — the single postulate and its consequences
- [Examples](examples/minkowski) — runnable Julia scripts
- [Documents](documents/) — 86 theoretical documents

---

*Working document. Speculative theoretical research.*
*© 2026 Jan Bouwman — MIT License*
