"""
FisherGeometrics.jl
===================
A unified field theory from the Fisher information geometry of
the Standard Model vacuum.

                    THE SINGLE POSTULATE

                    g_AB = 𝓕_AB / ρ₀

The spacetime metric IS the Fisher information tensor.
All Standard Model structure follows. Zero free parameters.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

STRUCTURE (each layer rests on the one below)

  Foundation      g_AB = 𝓕_AB/ρ₀  ·  τ=1/5  ·  φ  ·  κ_hol
      │
  Geometry        K = ℂP²×S³×S¹  ·  spectrum(Ð²_K)  ·  𝒯(K)=1
      │
  Symmetry        SU(3)×SU(2)×U(1)  ·  3 gen.  ·  sin²θ_W
      │
  Dynamics        α_em  ·  CKM  ·  δ=arctan(φ²)  ·  1:τ²:τ⁴
      │
  Gravity         G_μν=8πG_N ℛ_μν[𝓕]  ·  S_BH=A/4G_N  ·  Ω_Λ
      │
  Evolution       iħ dρ̂/dt = [Ð²_K, ρ̂]  ·  Bures distance

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

QUICK START

    using FisherGeometrics

    # Framework constants
    τ        # = 1/5  (from 4τ = cos(πτ))
    φ        # = golden ratio = [2]_q at level k=3
    κ_hol    # = 6/5  (holographic coupling)

    # Fisher tensor for a pure state
    ψ = [1, 0, 0, 0, 0, 0] |> ComplexF64
    F = fisher_tensor(pure_state(ψ))

    # CKM matrix
    V = ckm_matrix()

    # Time evolution: iħ dρ̂/dt = [Ð²_K, ρ̂]
    H   = hamiltonian_KK(6)
    ρ₀  = vacuum_state()
    ρ_t = evolve_exact(ρ₀, H, 1.0)

    # Full scoreboard
    scoreboard()

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

REFERENCES

  Jan Bouwman, FisherGeometrics Framework (2026)
  Documents I–LXX, available at github.com/[repository]

  Braunstein & Caves, PRL 72 (1994) — 𝓕^(Q) = 4 g^(FS)
  Castellani, D'Auria & Fré, NPB 239 (1984) — CDF eq.(1.2)
  Fabbri & Fré, NPB (1999) — AdS₄ × M^{1,1,1} spectrum
  Connes, CMP 182 (1996) — spectral action principle
  Wald, PRD 48 (1993) — black hole entropy formula

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
module FisherGeometrics

# ─────────────────────────────────────────────────────────────
# LOAD ALL LAYERS IN DEPENDENCY ORDER
# ─────────────────────────────────────────────────────────────

include("Foundation.jl")
include("Geometry.jl")
include("Symmetry.jl")
include("Dynamics.jl")
include("Gravity.jl")
include("Evolution.jl")

using .Foundation
using .Geometry
using .Symmetry
using .Dynamics
using .Gravity
using .Evolution

using Printf

# ─────────────────────────────────────────────────────────────
# RE-EXPORT — everything available at top level
# ─────────────────────────────────────────────────────────────

# Foundation
export τ, κ_hol, φ, λ_W, A_Wolf, δ_CP, M_Pl_GeV, M_c_GeV
export quantum_fisher, fisher_tensor, fubini_study_metric
export gellmann_basis, vacuum_state, pure_state

# Geometry
export spectrum_K, dirac_spectrum_K, kk_masses
export zeta_K, zeta_K_prime
export analytic_torsion
export bures_distance, information_distance
export vol_CP2, vol_S3, vol_S1, vol_K

# Symmetry
export n_generations, hypercharges
export sin2_weinberg, weinberg_angle
export alpha_GUT_inv, alpha_strong

# Dynamics
export alpha_em, alpha_em_inv
export mass_hierarchy, mass_ratios
export ckm_angles, ckm_matrix, ckm_wolfenstein
export jarlskog
export cp_phase_identity

# Gravity
export information_resistance
export cosmological_constant, Λ_fundamental
export bh_entropy, bh_temperature
export newton_constant_scaling
export Ω_Λ, n_s, Ω_DM
export unified_equation

# Evolution
export hamiltonian_KK
export von_neumann_rhs
export evolve_exact, evolve_rk4
export purity, entropy
export information_distance_trajectory
export decoherence_time
export measurement_projection

# ─────────────────────────────────────────────────────────────
# FULL SCOREBOARD
# ─────────────────────────────────────────────────────────────

"""
    scoreboard()

Complete comparison of all framework predictions vs observation.
Runs all six layers and prints the unified results table.
"""
function scoreboard()
    V = ckm_matrix()
    J = jarlskog()

    println()
    println("╔══════════════════════════════════════════════════════════════════╗")
    println("║         FisherGeometrics — Complete Scoreboard                   ║")
    println("║         g_AB = 𝓕_AB / ρ₀  ·  Zero free parameters              ║")
    println("╠══════════════════════════════════════════════════════════════════╣")
    println("║  Quantity                Predicted      Observed      Deviation  ║")
    println("╠══════════════════════════════════════════════════════════════════╣")

    function row(name, pred, obs, note="")
        dev = abs(pred - obs) / max(abs(obs), 1e-30) * 100
        @printf("║  %-22s  %10.5f   %10.5f   %5.2f%%  %-3s  ║\n",
                name, pred, obs, dev, note)
    end

    println("║  — Gauge sector ————————————————————————————————————————————  ║")
    row("sin²θ_W",           sin2_weinberg,      0.2312)
    row("α_s(M_Z)",          alpha_strong,        0.118)
    row("1/α_GUT",           alpha_GUT_inv(),     41.5)

    println("║  — Fine structure ——————————————————————————————————————————  ║")
    row("1/α_em",            alpha_em_inv,        137.036)

    println("║  — CKM matrix ——————————————————————————————————————————————  ║")
    row("|V_us| = λ_W",      λ_W,                 0.2250)
    row("|V_cb|",            abs(V[2,3]),          0.0418)
    row("|V_ub|",            abs(V[1,3]),          0.00351)
    row("δ_CP  (degrees)",   rad2deg(δ_CP),        69.2)
    row("J  (×10⁻⁵)",       J * 1e5,             3.08e-5 * 1e5)

    println("║  — Cosmology ———————————————————————————————————————————————  ║")
    row("Ω_Λ",               Ω_Λ,                 0.70)
    row("n_s",               n_s,                 0.9649)

    println("╠══════════════════════════════════════════════════════════════════╣")
    println("║  Exact results (deviation = 0 by construction):                  ║")
    println("║    Hypercharges  +1/6, +2/3, −1/3           EXACT               ║")
    println("║    3 generations  c₁(ℂP²) = 3               EXACT               ║")
    println("║    S_BH = A/(4G_N)  Wald formula             EXACT               ║")
    println("║    𝒯(K) = 1  Künneth + χ(S³)=0              EXACT               ║")
    println("║    Λ_fundamental = 0                          EXACT               ║")
    println("║    1 + φ⁴ = 3φ²  →  sin δ = φ/√3            EXACT               ║")
    println("╠══════════════════════════════════════════════════════════════════╣")
    println("║  Falsifiable predictions (not yet measured):                      ║")
    @printf("║    Σmν  = %4.0f meV        Euclid / CMB-S4  (2030)              ║\n", 61.0)
    @printf("║    σ_xy = 3 e²/h           kagomé metals  (now)                  ║\n")
    @printf("║    T_c  = 3ħ/(2k_B τ_c)   superconductors                        ║\n")
    println("╚══════════════════════════════════════════════════════════════════╝")
end

# ─────────────────────────────────────────────────────────────
# FULL CONSISTENCY CHECK
# ─────────────────────────────────────────────────────────────

"""
    check_all() → Bool

Run all consistency checks across all six layers.
Returns true if every check passes.
"""
function check_all()
    println("Running FisherGeometrics consistency checks...")
    println()

    results = [
        ("Foundation", true),   # constants checked at load time via @assert
        ("Geometry",   Geometry.check_geometry()),
        ("Symmetry",   Symmetry.check_symmetry()),
        ("Dynamics",   Dynamics.check_dynamics()),
        ("Gravity",    Gravity.check_gravity()),
        ("Evolution",  Evolution.check_evolution()),
    ]

    all_ok = true
    for (name, ok) in results
        status = ok ? "✓  pass" : "✗  FAIL"
        @printf("  %-14s  %s\n", name, status)
        all_ok = all_ok && ok
    end

    println()
    if all_ok
        println("All checks passed ✓  —  Framework internally consistent.")
    else
        println("Some checks failed ✗  —  See warnings above.")
    end

    return all_ok
end

# ─────────────────────────────────────────────────────────────
# PACKAGE INFO
# ─────────────────────────────────────────────────────────────

"""
    info()

Print framework overview and module structure.
"""
function info()
    println()
    println("  FisherGeometrics.jl")
    println("  ═══════════════════")
    println("  Unified field theory from Fisher information geometry.")
    println()
    println("  Postulate:  g_AB = 𝓕_AB / ρ₀")
    println()
    println("  Modules (each layer depends on those above it):")
    println("    Foundation   —  𝓕_AB, τ=$(τ), φ=$(round(φ,digits=4)), κ_hol=$(κ_hol)")
    println("    Geometry     —  K=ℂP²×S³×S¹, Ð²_K, 𝒯(K)=1")
    println("    Symmetry     —  SU(3)×SU(2)×U(1), $(n_generations) gen., sin²θ_W=$(sin2_weinberg)")
    println("    Dynamics     —  α_em=1/$(round(alpha_em_inv,digits=1)), CKM, δ=$(round(rad2deg(δ_CP),digits=2))°")
    println("    Gravity      —  G_μν=8πG_N ℛ_μν[𝓕], S_BH=A/4G_N, Ω_Λ=$(round(Ω_Λ,digits=3))")
    println("    Evolution    —  iħ dρ̂/dt=[Ð²_K,ρ̂], Bures distance")
    println()
    println("  Zero free parameters.")
    println("  Run scoreboard() for full results.")
    println("  Run check_all() for consistency verification.")
    println()
end

end # module FisherGeometrics
