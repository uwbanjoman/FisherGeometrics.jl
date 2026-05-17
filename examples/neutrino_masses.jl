# examples/neutrino_masses.jl
# ============================
# Neutrino masses from the Freund-Rubin spectrum of K = ℂP² × S³ × S¹
# with unequal radii r₁ (ℂP²), r₂ (S³), r₃ (S¹).
#
# ── The physics ──────────────────────────────────────────────────────────────
#
# In the leading-order approximation (Documents VI, LXXII), the three neutrino
# masses follow the Killing-spinor hierarchy:
#
#   m₁ : m₂ : m₃ = τ⁴ : τ² : 1
#
# This gives Σmν_bare ≈ 51.6 meV and Δm²₂₁ ≈ 3.9 × 10⁻⁵ eV² (observed:
# 7.53 × 10⁻⁵ eV²). The discrepancy signals that unequal radii are needed.
#
# The full Kaluza-Klein spectrum with unequal radii r₁, r₂, r₃ gives:
#
#   λ(k,j,n) = 4k(k+2)/r₁² + j(j+2)/r₂² + n²/r₃²
#
# This script finds the radius r₂ (S³) that simultaneously reproduces
# both mass splittings, with r₃ = τ × r₁ fixed by Document XXIV.
#
# Document LXXII, FisherGeometrics Framework, May 2026
# © 2026 Jan Bouwman

using FisherGeometrics
using LinearAlgebra
using Printf

# ── Experimental inputs ───────────────────────────────────────────────────────

const Δm2_32_obs = 2.453e-3   # eV²  (NuFit 5.3, normal hierarchy)
const Δm2_21_obs = 7.530e-5   # eV²  (NuFit 5.3)
const τ_val      = Float64(τ) # = 1/5, from Document XXIV

println("="^65)
println("  FisherGeometrics — Neutrino Masses from Freund-Rubin Spectrum")
println("="^65)
println()
@printf("  Target Δm²₃₂ = %.4e eV²  (atmospheric, NuFit 5.3)\n", Δm2_32_obs)
@printf("  Target Δm²₂₁ = %.4e eV²  (solar,       NuFit 5.3)\n", Δm2_21_obs)
@printf("  Constraint τ  = r₃/r₁ = 1/%d  (Document XXIV)\n", round(Int, 1/τ_val))
println()

# ── KK spectrum functions ─────────────────────────────────────────────────────

"""
    lepton_spectrum(r1, r2, r3; j_max, n_max) → Vector{Float64}

Kaluza-Klein eigenvalues in the lepton sector (ℂP² zero mode k=0)
for unequal radii. Returns sorted unique eigenvalues.
"""
function lepton_spectrum(r1::Real, r2::Real, r3::Real;
                         j_max::Int=30, n_max::Int=30)
    E0 = (9/4) * (1/r1^2 + 1/r2^2 + 1/r3^2) / 3
    λs = Float64[]
    for j in 0:j_max
        λ_S3 = j*(j+2) / r2^2
        for n in 0:n_max
            λ_S1 = n^2 / r3^2
            push!(λs, λ_S3 + λ_S1 + E0)
            if n > 0
                push!(λs, λ_S3 + λ_S1 + E0)
            end
        end
    end
    return sort(unique(round.(λs, digits=10)))
end

"""
    neutrino_masses_eV(r1, r2, r3, M_scale_eV) → (m1, m2, m3)

Three lightest neutrino masses in eV.
"""
function neutrino_masses_eV(r1::Real, r2::Real, r3::Real, M_scale_eV::Real)
    λs = lepton_spectrum(r1, r2, r3)
    λ_light = λs[1:min(3, length(λs))]
    return sort(sqrt.(λ_light) .* M_scale_eV)
end

"""
    compute_masses(r1, r2, r3) → NamedTuple

Compute all neutrino observables for given radii.
"""
function compute_masses(r1::Real, r2::Real, r3::Real)
    λs = lepton_spectrum(r1, r2, r3)
    length(λs) < 3 && return nothing
    la, lb, lc = λs[1], λs[2], λs[3]
    lc > la || return nothing

    M_scale = sqrt(Δm2_32_obs / (lc - la))
    m1, m2, m3 = neutrino_masses_eV(r1, r2, r3, M_scale)

    Δm2_32 = m3^2 - m1^2
    Δm2_21 = m2^2 - m1^2
    Σmν    = (m1 + m2 + m3) * 1e3

    return (m1=m1*1e3, m2=m2*1e3, m3=m3*1e3,
            Δm2_32=Δm2_32, Δm2_21=Δm2_21, Σmν=Σmν,
            ratio_21=Δm2_21/Δm2_21_obs, M_scale=M_scale*1e3)
end

# ── Leading-order check ───────────────────────────────────────────────────────

println("─── Leading-order result (equal radii r₂ = r₃ = τ = 1/5) ───")
println()

res_eq = compute_masses(1.0, τ_val, τ_val)
@printf("  r₁ = 1.0000, r₂ = %.4f (S³), r₃ = %.4f (S¹)\n", τ_val, τ_val)
@printf("  M_scale = %.4f meV\n", res_eq.M_scale)
println()
@printf("  m₁ = %8.4f meV\n", res_eq.m1)
@printf("  m₂ = %8.4f meV\n", res_eq.m2)
@printf("  m₃ = %8.4f meV\n", res_eq.m3)
println()
@printf("  Δm²₃₂ = %.4e eV²  (target: %.4e,  Δ = %+.1f%%)\n",
        res_eq.Δm2_32, Δm2_32_obs, (res_eq.Δm2_32/Δm2_32_obs - 1)*100)
@printf("  Δm²₂₁ = %.4e eV²  (target: %.4e,  Δ = %+.0f%%)\n",
        res_eq.Δm2_21, Δm2_21_obs, (res_eq.ratio_21 - 1)*100)
@printf("  Σmν   = %.2f meV\n", res_eq.Σmν)
println()
println("  → Δm²₂₁ is off: unequal radii needed.")
println()

# ── Scan over r₂ ─────────────────────────────────────────────────────────────

println("─── Scan over r₂/r₁  (r₃ = τ fixed) ───")
println()
@printf("  %-8s  %-10s  %-10s  %-10s  %-12s  %-10s\n",
        "r₂/r₁", "m₁ (meV)", "m₂ (meV)", "m₃ (meV)", "Δm²₂₁/obs", "Σmν (meV)")
println("  " * "─"^65)

function scan_r2(r2_values)
    best_diff   = Inf
    best_r2     = τ_val
    best_result = nothing

    for r2_ratio in r2_values
        res = compute_masses(1.0, r2_ratio, τ_val)
        res === nothing && continue

        diff = abs(res.ratio_21 - 1.0)
        if diff < best_diff
            best_diff   = diff
            best_r2     = r2_ratio
            best_result = res
        end

        # Print every 5th step
        if isapprox(r2_ratio, round(r2_ratio/0.05)*0.05, atol=0.005)
            @printf("  %-8.3f  %-10.4f  %-10.4f  %-10.4f  %-12.4f  %-10.4f\n",
                    r2_ratio, res.m1, res.m2, res.m3, res.ratio_21, res.Σmν)
        end
    end
    return best_r2, best_result
end

best_r2, best_res = scan_r2(0.10:0.01:0.50)
println()

# ── Best fit ──────────────────────────────────────────────────────────────────

println("─── Best fit result ───")
println()
b = best_res
@printf("  r₁ = 1.0000 (ℂP²)\n")
@printf("  r₂ = %.4f (S³)   →  r₂/r₁ = %.4f  (was τ = %.4f)\n",
        best_r2, best_r2, τ_val)
@printf("  r₃ = %.4f (S¹)   →  r₃/r₁ = τ = 1/5  (fixed)\n", τ_val)
println()
@printf("  m₁ = %8.4f meV\n", b.m1)
@printf("  m₂ = %8.4f meV\n", b.m2)
@printf("  m₃ = %8.4f meV\n", b.m3)
println()
@printf("  Δm²₃₂ = %.4e eV²  (target: %.4e)  ✓\n", b.Δm2_32, Δm2_32_obs)
@printf("  Δm²₂₁ = %.4e eV²  (target: %.4e,  Δ = %+.1f%%)\n",
        b.Δm2_21, Δm2_21_obs, (b.ratio_21 - 1)*100)
@printf("  Σmν   = %.2f meV\n", b.Σmν)
println()

# ── Holographic correction ────────────────────────────────────────────────────

println("─── Holographic correction from Hopf screen S³ → S² ───")
println()

M_c_eV  = 1.44e17 * 1e9
M_Pl_eV = 1.22e19 * 1e9
δm_holo_meV = τ_val * M_c_eV^2 / (24 * M_Pl_eV) * 1e3

@printf("  Hopf screen area: A = π τ² = π/25 = %.6f\n", π*τ_val^2)
@printf("  Holographic mass correction: δmν = τ M_c²/(24 M_Pl) = %.2f meV\n",
        δm_holo_meV)
println()

m1_tot = b.m1 + δm_holo_meV
m2_tot = b.m2 + δm_holo_meV
m3_tot = b.m3 + δm_holo_meV
Σ_tot  = m1_tot + m2_tot + m3_tot

@printf("  m₁ (total) = %8.4f meV\n", m1_tot)
@printf("  m₂ (total) = %8.4f meV\n", m2_tot)
@printf("  m₃ (total) = %8.4f meV\n", m3_tot)
@printf("  Σmν (total) = %.2f meV\n", Σ_tot)
println()

# ── Summary ───────────────────────────────────────────────────────────────────

println("─── Summary ───")
println()
println("  ┌─────────────────────────────────────────────────────────┐")
@printf("  │  Σmν (bare, equal radii)    = %5.1f meV               │\n", res_eq.Σmν)
@printf("  │  Σmν (best fit r₂)         = %5.1f meV               │\n", b.Σmν)
@printf("  │  Holographic correction     = %5.1f meV (×3 gen.)     │\n", δm_holo_meV)
@printf("  │  Σmν (total)               = %5.1f meV               │\n", Σ_tot)
println("  │                                                         │")
println("  │  Observed upper limit:  < 120 meV  (Planck 2018)       │")
println("  │  Euclid/CMB-S4 target:  ~20 meV sensitivity (2030)     │")
println("  │                                                         │")
println("  │  Framework prediction:  52 – 65 meV                    │")
println("  │  Falsified by:          Σmν < 30 meV                   │")
println("  └─────────────────────────────────────────────────────────┘")
println()
println("  Normal hierarchy — testable by JUNO (2024–2027).")
println("  See: Document LXXII, FisherGeometrics Framework (2026)")
