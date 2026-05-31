---
title: Internal Space
parent: Framework
nav_order: 2
---

# The Internal Space K = ℂP² × S³ × S¹

## Why K is not a postulate

K is derived, not assumed. The derivation uses three steps
(Document LXXXVI):

**Step 1 — Minimality**

The minimal Hilbert space carrying colour, isospin, and hypercharge
simultaneously is ℋ = ℂ⁶ = ℂ³ ⊗ ℂ².

This follows from the dimension formula for irreducible representations:
- SU(3) needs at least ℂ³ (no rep with 1 < dim < 3 exists)
- SU(2) needs at least ℂ² (dim = 2j+1, minimum non-trivial j=½)

**Step 2 — Anomaly freedom**

The unique hypercharge assignment making 𝓕_AB gauge-invariant on ℂ⁶:

$$\sum Y^3 = 0, \quad \sum Y = 0$$

gives Y = (+1/6, +2/3, −1/3, −1/2, −1). Verified symbolically
with HomotopyContinuation.jl (`examples/derive_K.jl`).

**Step 3 — Sector factorisation**

The Fisher tensor factorises into three independent sectors:

$$K = K_{\rm colour} \times K_{\rm weak} \times K_Y
  = \mathbb{CP}^2 \times S^3 \times S^1$$

---

## The parameter τ = 1/5

The single geometric parameter of K — the ratio of radii
τ = r_{S¹}/r_{ℂP²} — is fixed by two independent equations:

$$\kappa_{\rm hol} = 6\tau \qquad \text{(information density)}$$

$$\kappa_{\rm hol} = \frac{3}{2}(1-\tau) \qquad \text{(Killing spinor integrability)}$$

Setting equal: **τ = 1/5**, κ_hol = 6/5.

The golden ratio follows: φ = 2cos(πτ) = (1+√5)/2.

---

## The honest epistemological status

This has the same status as general relativity:

> GR does not explain why there are 3+1 spacetime dimensions,
> but *given* 3+1 dimensions and the equivalence principle,
> the Einstein equation is uniquely determined.

Here: the framework does not explain why nature has colour, isospin,
and hypercharge. But *given* these three quantum numbers,
K is uniquely determined.

---

## Symbolic derivation

The full derivation runs in seconds:

```bash
julia --project examples/derive_K.jl
```

Uses Symbolics.jl + HomotopyContinuation.jl to solve the anomaly
equations and Fisher self-consistency equations symbolically.
