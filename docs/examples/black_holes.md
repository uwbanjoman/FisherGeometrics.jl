---
title: Black Holes
parent: Examples
nav_order: 3
---

# Black Holes and the Information Paradox

Demonstrates that black holes have no singularity and no information
paradox in the FisherGeometrics framework.

## Run it

```bash
julia --project examples/black_holes.jl
```

## What it demonstrates

**No singularity**

In classical GR, the metric diverges at r=0. In FisherGeometrics,
the fundamental object ρ̂ is always a bounded 6×6 matrix.
The "singularity" is an artefact of using the metric as the
fundamental object.

```
Black hole centre:  λ = [1.000, 0.000, 0.000, ...]  Tr(ρ²)=1.000 ✓
Adjacent site:      λ = [0.750, 0.050, 0.050, ...]  Tr(ρ²)=0.575 ✓
Far from BH:        λ = [0.167, 0.167, 0.167, ...]  Tr(ρ²)=0.167 ✓
```

**No information paradox**

The Von Neumann equation is unitary. Information is never destroyed.
The "paradox" does not arise when ρ̂ is the fundamental object.

**Hawking evaporation**

$$\hat\rho_{\rm BH}(t) \to \mathbf{I}/6 \quad (t\to\infty)$$

The black hole state evolves from a pure state toward the vacuum.
Entropy increases monotonically — the Page curve is automatic.

## The event horizon

The horizon is where the informative velocity vanishes:

$$v_{\rm info} = \left|\frac{d\hat\rho}{dt}\right|_{\mathcal{F}} \to 0$$

Information flow freezes at the horizon. But ρ̂ remains smooth.

## Related documents

- Document I — Bekenstein-Hawking entropy S_BH = A/4G_N
- Phase 5 — Black hole simulation
