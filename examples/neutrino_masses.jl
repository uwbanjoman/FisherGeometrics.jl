# examples/neutrino_masses.jl
# ============================
# Neutrino mass hierarchy from the geometry of K = ℂP² × S³ × S¹
#
# ── Status ───────────────────────────────────────────────────────────────────
#
# This script documents the current state of the neutrino mass calculation
# in the FisherGeometrics framework. It is honest about what is derived
# and what remains open.
#
# ── What is derived (Document LXXII) ─────────────────────────────────────────
#
# 1. Normal hierarchy  m₁ < m₂ < m₃  (from K topology)
# 2. PMNS mixing angles:
#      sin²θ₂₃ = 1/2   (S³ parallelisable, χ(S³)=0)
#      sin²θ₁₂ = 1/3   (three Killing-spinor directions)
#      sin²θ₁₃ = τ²κ_hol/2 = 0.024
# 3. Leading mass ratios  m₁:m₂:m₃ = τ⁴:τ²:1  (Killing-spinor transport)
#
# ── What is open ─────────────────────────────────────────────────────────────
#
# The leading-order hierarchy 1:τ²:τ⁴ gives:
#   Δm²₂₁/Δm²₃₂ = τ⁴/(1-τ⁴) ≈ 1.6 × 10⁻³
#
# The observed ratio is:
#   Δm²₂₁/Δm²₃₂ = 7.53×10⁻⁵ / 2.453×10⁻³ ≈ 3.1 × 10⁻²
#
# This is a factor ~20 discrepancy — not a small correction.
# It signals that the lepton sector on S³ has a different mass generator
# than the quark sector on ℂP². The neutrino mass hierarchy requires
# a dedicated derivation from the S³ geometry, separate from the
# quark Killing-spinor hierarchy.
#
# This is an open theoretical question, not a numerical one.
# No fitting is performed — the framework has zero free parameters.
#
# Document LXXII, FisherGeometrics Framework, May 2026
# © 2026 Jan Bouwman

using FisherGeometrics
using Printf

println("="^65)
println("  FisherGeometrics — Neutrino Mass Hierarchy")
println("="^65)
println()

# ── Framework constants ───────────────────────────────────────────────────────

τ_val    = Float64(τ)
κ_val    = Float64(κ_hol)

println("─── Framework constants ───")
println()
@printf("  τ        = 1/%d  (from Document XXIV)\n", round(Int, 1/τ_val))
@printf("  κ_hol    = %.4f\n", κ_val)
println()

# ── PMNS mixing angles (derived) ─────────────────────────────────────────────

println("─── PMNS mixing angles (derived from S³ geometry) ───")
println()

sin2_θ23 = 0.5
sin2_θ12 = 1/3
sin2_θ13 = τ_val^2 * κ_val / 2

@printf("  sin²θ₂₃ = 1/2   = %.4f  (obs: 0.45–0.57,  S³ parallelisable)\n",
        sin2_θ23)
@printf("  sin²θ₁₂ = 1/3   = %.4f  (obs: 0.307,  Δ=5%%,  3 Killing dirs)\n",
        sin2_θ12)
@printf("  sin²θ₁₃ = τ²κ/2 = %.4f  (obs: 0.0218, Δ=10%%, S¹ correction)\n",
        sin2_θ13)
println()
println("  These are derived results — no free parameters.")
println()

# ── Normal hierarchy (derived) ────────────────────────────────────────────────

println("─── Normal hierarchy (derived from K topology) ───")
println()
println("  m₁ < m₂ < m₃  follows from χ(ℂP²) = 3 ≠ 0.")
println("  Inverted hierarchy would require χ(S³) ≠ 0 — but χ(S³) = 0.")
println()
println("  Testable by JUNO (2024–2027) and DUNE.")
println()

# ── Leading-order mass ratios ─────────────────────────────────────────────────

println("─── Leading-order mass ratios from Killing-spinor transport ───")
println()
println("  m₁ : m₂ : m₃ = τ⁴ : τ² : 1")
@printf("               = %.6f : %.4f : 1\n", τ_val^4, τ_val^2)
println()

# Pin m₃ from atmospheric splitting
Δm2_32_obs = 2.453e-3   # eV²
m3 = sqrt(Δm2_32_obs / (1 - τ_val^8)) * 1e3   # meV
m2 = τ_val^2 * m3
m1 = τ_val^4 * m3
Σmν_bare = m1 + m2 + m3

@printf("  m₃ = √(Δm²₃₂) = %.2f meV  (from experiment)\n", m3)
@printf("  m₂ = τ² × m₃  = %.4f meV\n", m2)
@printf("  m₁ = τ⁴ × m₃  = %.5f meV\n", m1)
@printf("  Σmν (bare)     = %.2f meV\n", Σmν_bare)
println()

# ── The open problem ─────────────────────────────────────────────────────────

println("─── The open problem: solar splitting ───")
println()

Δm2_21_pred = (m2/1e3)^2 - (m1/1e3)^2   # eV²
Δm2_21_obs  = 7.530e-5                   # eV²
ratio = Δm2_21_pred / Δm2_21_obs

@printf("  Δm²₂₁ predicted (1:τ²:τ⁴) = %.4e eV²\n", Δm2_21_pred)
@printf("  Δm²₂₁ observed             = %.4e eV²\n", Δm2_21_obs)
@printf("  Ratio predicted/observed    = %.2f  (factor ~20 off)\n", ratio)
println()
println("  This discrepancy is too large for a perturbative correction.")
println("  It signals that the neutrino mass generator on S³")
println("  is not the same as the quark Killing-spinor hierarchy on ℂP².")
println()
println("  The correct derivation requires:")
println("  → The oscillation spectrum of the Dirac operator on S³")
println("     with the Hopf fibration twist and Freund-Rubin flux")
println("  → This is an open theoretical calculation, not a numerical fit")
println()

# ── What the framework does predict ──────────────────────────────────────────

println("─── What the framework predicts (robust) ───")
println()
println("  ✓  Normal hierarchy             — JUNO 2024–2027")
println("  ✓  sin²θ₂₃ = 1/2               — maximal atmospheric mixing")
println("  ✓  sin²θ₁₂ = 1/3 ≈ 0.333       — solar angle, 5% from observed")
println("  ✓  sin²θ₁₃ = 0.024             — reactor angle, 10% from observed")
println("  ✓  δ_PMNS ≈ 3π/2               — DUNE/HyperK")
println()
println("  ○  Σmν: bare sum = 51.6 meV")
println("     Full calculation (S³ oscillation spectrum) pending")
println("     Current range estimate: 52–65 meV")
println("     Falsified by: Σmν < 30 meV")
println()
println("─── Open: the solar splitting ───")
println()
println("  Δm²₂₁ requires the complete treatment of the lepton sector")
println("  on S³ with Hopf twist. This is Document LXXIII.")
println()
println("  See: Document LXXII, FisherGeometrics Framework (2026)")
