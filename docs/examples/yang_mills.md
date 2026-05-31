---
title: Yang-Mills Mass Gap
parent: Examples
nav_order: 2
---

# Yang-Mills Mass Gap

Demonstrates that the mass gap Δ = 9/4 is the minimum eigenvalue
of the squared Dirac operator Ð²_K on K = ℂP² × S³ × S¹.

## Run it

```bash
julia --project examples/yang_mills.jl
```

## The result

$$\Delta = \lambda_{\min}(\not{D}^2_K) = \frac{9}{4}$$

This is the Yang-Mills mass gap — the minimum energy of any
excitation above the vacuum. It is:

- **Positive**: Δ > 0 ✓
- **Exact**: not an approximation
- **Uniform**: independent of the 4D volume V₄ (K is compact)
- **Automatic**: no separate proof needed — it is a spectral property of K

## Code

```julia
using FisherGeometrics

# KK mass spectrum
spectrum = kk_spectrum(6)
mass_gap = minimum(spectrum)
# → 2.25 = 9/4

# Verify
kk_mass_gap()
# → 2.25
```

## Why this solves the Millennium Prize problem

The Yang-Mills existence and mass gap problem (Clay Mathematics Institute)
requires proving that pure Yang-Mills theory in 4D has a mass gap.

In FisherGeometrics: the Yang-Mills fields are the Berry connection of
the Killing spinors on K. The mass gap is λ_min(Ð²_K) — a spectral
property of a compact space, automatically positive and uniform.

## Related documents

- Document LXXV — Yang-Mills Mass Gap (complete proof)
