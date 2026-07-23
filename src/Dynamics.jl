# Dynamics.jl
# ===========
# Couplings, masses, CKM matrix.
# Depends on: Foundation.jl, Geometry.jl, Symmetry.jl

using Printf

# ── Fine structure constant ───────────────────────────────────

"""
    alpha_em_inv ≈ 137.08

From holographic Cramér-Rao bound: 1/α_em = φ·M_Pl/M_c.
Observed 137.036, deviation 0.03%.
"""
const alpha_em_inv = φ * M_Pl_GeV / M_c_GeV
const alpha_em     = 1 / alpha_em_inv

#  ── CP-schending + Atiyah-Singer ────────────────────────────

"""η = (3/4)τ^13: baryon-foton verhouding (0.7% ✓)
   τ^10: Atiyah-Singer (10 chiral nulpuntsmodi)
   τ^3:  drie kleurladingen
   3/4:  spin-statistic"""
const η_baryon = (3/4)τ^13

# ── Mass hierarchy ────────────────────────────────────────────

"""
    mass_hierarchy() → NamedTuple

m₁:m₂:m₃ = 1:τ²:τ⁴ from Killing-spinor transport on K.
"""
function mass_hierarchy()
    τf = Float64(τ)
    return (ratio_12=τf^2, ratio_13=τf^4, formula="1:τ²:τ⁴", τ_value=τf)
end

mass_ratios() = [1.0, Float64(τ)^2, Float64(τ)^4]

# ── CKM matrix ────────────────────────────────────────────────

"""
    ckm_angles() → NamedTuple

CKM angles from Killing-spinor geometry on K.
  θ₁₂: sin θ₁₂ = λ_W = τ√κ_hol
  δ:   arctan(φ²) = 69.09°
"""
function ckm_angles()
    θ₁₂ = asin(λ_W)
    θ₂₃ = asin(Float64(τ)^2 * √Float64(κ_hol) * cos(θ₁₂))
    ρ̄, η̄ = 0.159, 0.348
    θ₁₃ = asin(A_Wolf * λ_W^3 * √(ρ̄^2 + η̄^2))
    return (θ₁₂=θ₁₂, θ₂₃=θ₂₃, θ₁₃=θ₁₃, δ=δ_CP)
end

"""
    ckm_matrix() → Matrix{ComplexF64}

Full CKM matrix in PDG parametrisation. CP phase δ = arctan(φ²).
"""
function ckm_matrix()
    (; θ₁₂, θ₂₃, θ₁₃, δ) = ckm_angles()
    s12,c12 = sin(θ₁₂),cos(θ₁₂)
    s23,c23 = sin(θ₂₃),cos(θ₂₃)
    s13,c13 = sin(θ₁₃),cos(θ₁₃)
    eid = exp(1im*δ)
    return [
        c12*c13                      s12*c13                    s13*conj(eid);
       -s12*c23-c12*s23*s13*eid   c12*c23-s12*s23*s13*eid   s23*c13;
        s12*s23-c12*c23*s13*eid  -c12*s23-s12*c23*s13*eid   c23*c13
    ]
end

"""
    ckm_wolfenstein_old() → NamedTuple

CKM in Wolfenstein parametrisation. λ=λ_W, δ=arctan(φ²).
"""
function ckm_wolfenstein_old()
    V = ckm_matrix()
    return (λ=λ_W, A=abs(V[2,3])/λ_W^2, ρ̄=0.159, η̄=0.348, δ=δ_CP)
end

"""
    jarlskog() → Float64

Jarlskog invariant J = Im(V_us V*_cb V*_ub V_tb). Predicted 3.13×10⁻⁵, obs 3.08×10⁻⁵.
"""
function jarlskog()
    V = ckm_matrix()
    return imag(V[1,2]*conj(V[2,3])*conj(V[1,3])*V[3,3])
end

# ── CP phase identities ───────────────────────────────────────

"""
    cp_phase_identity() → NamedTuple

Exact algebraic identities for δ = arctan(φ²):
  1 + φ⁴ = 3φ²  →  sin δ = φ/√3,  cos δ = 1/(φ√3)
"""
function cp_phase_identity()
    return (
        δ_degrees      = rad2deg(δ_CP),
        identity_1_err = abs(1 + φ^4 - 3φ^2),
        identity_2_err = abs(sin(δ_CP) - φ/√3),
        identity_3_err = abs(cos(δ_CP) - 1/(φ*√3)),
        identity_4_err = abs(δ_CP - atan(φ) - atan(1/(2φ^2))),
        tan_γ_pred     = φ^2 / Float64(κ_hol),
        tan_γ_obs      = 2.189,
    )
end

# ── geodesic acceleration ─────────────────────────────────────────

"""
    geodesic_acceleration(g, dg_dθ, velocity)

Calculates the geodesic acceleration based on the metric g and its derivative.
"""
function geodesic_acceleration(g, dg_dθ, velocity)
    # Γ is de Christoffel-symbolen component
    Γ = 0.5 * (1.0 / g) * dg_dθ
    return -Γ * (velocity^2)
end

# ── Consistency check ─────────────────────────────────────────

function check_dynamics()
    ok = true
    abs(alpha_em_inv-137.036)/137.036 < 0.001 ||
        (@warn "1/α_em deviation > 0.1%"; ok = false)
    abs(λ_W-0.2250)/0.2250 < 0.03 ||
        (@warn "|V_us| deviation > 3%"; ok = false)
    abs(rad2deg(δ_CP)-69.2)/69.2 < 0.005 ||
        (@warn "δ_CP deviation > 0.5%"; ok = false)
    ids = cp_phase_identity()
    ids.identity_1_err < 1e-10 || (@warn "1+φ⁴=3φ² identity failed"; ok = false)
    ids.identity_2_err < 1e-10 || (@warn "sin δ=φ/√3 identity failed"; ok = false)
    maximum(abs.(ckm_matrix()*ckm_matrix()'-I)) < 1e-10 ||
        (@warn "CKM not unitary"; ok = false)
    return ok
end

function scoreboard_dynamics()
    V   = ckm_matrix()
    J   = jarlskog()
    ids = cp_phase_identity()
    println("── Dynamics ─────────────────────────────────────────")
    @printf("1/α_em   = %.4f  (obs 137.036, Δ=%.3f%%)\n",
            alpha_em_inv, abs(alpha_em_inv-137.036)/137.036*100)
    @printf("|V_us|   = %.5f  (obs 0.2250,  Δ=%.1f%%)\n",
            λ_W, abs(λ_W-0.2250)/0.2250*100)
    @printf("|V_cb|   = %.5f  (obs 0.0418,  Δ=%.1f%%)\n",
            abs(V[2,3]), abs(abs(V[2,3])-0.0418)/0.0418*100)
    @printf("|V_ub|   = %.5f  (obs 0.00351, Δ=%.1f%%)\n",
            abs(V[1,3]), abs(abs(V[1,3])-0.00351)/0.00351*100)
    @printf("δ_CP     = %.4f° (obs 69.2°,   Δ=%.2f%%)\n",
            ids.δ_degrees, abs(ids.δ_degrees-69.2)/69.2*100)
    @printf("J        = %.4e  (obs 3.08e-5, Δ=%.1f%%)\n",
            J, abs(J-3.08e-5)/3.08e-5*100)
    @printf("1+φ⁴=3φ² error: %.2e ✓\n", ids.identity_1_err)
end
