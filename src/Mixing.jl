# src/Mixing.jl

"""
    ckm_wolfenstein(; verbose::Bool=true) -> NamedTuple

CKM quark mixing matrix in the Wolfenstein parametrisation,
derived from the geometry of CP² × S³.

All four parameters follow from τ = 1/5 and the golden ratio
φ = (1+√5)/2 — zero free parameters:

    λ     = 9τ/8          SU(3)/SU(2) Clebsch-Gordan structure
    A     = 4τ + τ²       S³ spin-3/2 coupling to b-quark
    η̄     = 7τ/4          Hopf fibration S³ → S²
    δ_CP  = arctan(φ²)    Chern-Simons structure of K (universal)

The Cabibbo parameter λ = 9τ/8 = 0.22500 is exact to five decimal
places. The CP-violating phase δ_CP = arctan(φ²) = 69.09° is
identical for quarks and neutrinos — both determined by the same
Hopf fibration of the internal space K.

Average deviation from PDG 2024: 0.64%.

# Arguments
- `verbose`: print the parameter table (default: true)

# Returns
NamedTuple with fields: λ, A, η̄, δ_CP, J (Jarlskog invariant)

# Examples
```julia
ckm = ckm_wolfenstein()
ckm.λ       # → 0.22500  (exact)
ckm.δ_CP    # → 69.095°
ckm.J       # → Jarlskog invariant
```

See: FisherGeometrics preprint v15, section 8.
     examples/ckm_derivation.jl for the full CKM matrix.
"""
function ckm_wolfenstein(; verbose::Bool=true)
    τ = 1/5
    φ = (1+sqrt(5))/2

    λ_FG   = 9τ/8
    A_FG   = 4τ + τ^2
    η̄_FG   = 7τ/4
    δ_FG   = rad2deg(atan(φ^2))
    J_FG   = A_FG^2 * λ_FG^6 * η̄_FG

    if verbose
        λ_obs = 0.22500; A_obs = 0.826
        η̄_obs = 0.348;   δ_obs = 69.3

        println("\nCKM WOLFENSTEIN PARAMETERS — FisherGeometrics")
        println("="^60)
        @printf("%-8s  %-10s  %-10s  %-8s  %-20s\n",
                "Param", "FG", "PDG 2024", "Error", "Formula")
        println("─"^60)
        for (name, pred, obs, formula) in [
            ("λ",     λ_FG,  λ_obs, "9τ/8"),
            ("A",     A_FG,  A_obs, "4τ + τ²"),
            ("η̄",     η̄_FG,  η̄_obs, "7τ/4"),
            ("δ_CP°", δ_FG,  δ_obs, "arctan(φ²)"),
        ]
            @printf("%-8s  %-10.5f  %-10.5f  %-8.2f%%  %s\n",
                    name, pred, obs, abs(pred-obs)/obs*100, formula)
        end
        println()
        @printf("  Jarlskog invariant J = %.4e\n", J_FG)
        @printf("  Note: δ_CP (quarks) = δ_CP (neutrinos) = arctan(φ²) ✓\n")
        println()
    end

    return (λ=λ_FG, A=A_FG, η̄=η̄_FG, δ_CP=δ_FG, J=J_FG)
end

"""
    pmns_angles(; verbose::Bool=true) -> NamedTuple

PMNS neutrino mixing matrix angles from the geometry of CP² × S³ × S¹.

All four parameters follow from τ = 1/5 and φ = (1+√5)/2 — zero free
parameters:

    θ₁₂ = arctan(1/√2)      geodesic distance in CP¹ ⊂ CP²
                              (Fubini-Study metric, tribimaximal limit)
    θ₁₃ = arcsin(3τ/4)      S¹ fiber radius τ relative to CP²
    θ₂₃ = π/4 + arctan(τ²)  SU(2) maximal mixing + τ² correction
    δ_CP = arctan(φ²)        Hopf fibration S³ → S² (universal)

The CP-violating phase δ_CP = arctan(φ²) = 69.09° is identical for
quarks (CKM) and neutrinos (PMNS) — both determined by the same
Hopf fibration of the internal space K = CP² × S³ × S¹.

CP violation is not a free parameter but a topological invariant of K.

Average deviation from PDG 2024: 2.5%.
Unitarity: ‖UU† − I‖ = 1.85×10⁻¹⁶ ✓

# Arguments
- `verbose`: print the angle table (default: true)

# Returns
NamedTuple with fields: θ₁₂, θ₁₃, θ₂₃, δ_CP (all in degrees)

# Examples
```julia
pmns = pmns_angles()
pmns.θ₁₃     # → 8.63°   (measured: 8.57°, 0.7%)
pmns.δ_CP    # → 69.09°  (measured: 69.2°, 0.2%)
```

See: FisherGeometrics preprint v15, section 9.
     examples/pmns_derivation.jl for the full PMNS matrix.
"""
function pmns_angles(; verbose::Bool=true)
    τ = 1/5
    φ = (1+sqrt(5))/2

    θ₁₂ = rad2deg(atan(1/sqrt(2)))
    θ₁₃ = rad2deg(asin(3τ/4))
    θ₂₃ = rad2deg(π/4 + atan(τ^2))
    δ_CP = rad2deg(atan(φ^2))

    if verbose
        println("\nPMNS MIXING ANGLES — FisherGeometrics")
        println("="^60)
        @printf("%-8s  %-10s  %-10s  %-8s  %-22s\n",
                "Angle", "FG", "PDG 2024", "Error", "Formula")
        println("─"^60)
        for (name, pred, obs, formula) in [
            ("θ₁₂°", θ₁₂, 33.44, "arctan(1/√2)"),
            ("θ₁₃°", θ₁₃,  8.57, "arcsin(3τ/4)"),
            ("θ₂₃°", θ₂₃, 49.20, "π/4 + arctan(τ²)"),
            ("δ_CP°", δ_CP, 69.20, "arctan(φ²)"),
        ]
            @printf("%-8s  %-10.4f  %-10.4f  %-8.2f%%  %s\n",
                    name, pred, obs, abs(pred-obs)/obs*100, formula)
        end
        println()
        @printf("  Note: δ_CP (quarks) = δ_CP (neutrinos) = arctan(φ²) ✓\n")
        @printf("  Both determined by the Hopf fibration S³ → S²\n")
        println()
    end

    return (θ₁₂=θ₁₂, θ₁₃=θ₁₃, θ₂₃=θ₂₃, δ_CP=δ_CP)
end

"""
    cp_phase() -> Float64

The universal CP-violating phase from the Hopf fibration S³ → S².

    δ_CP = arctan(φ²)  where φ = (1+√5)/2 (golden ratio)

This phase is identical for quarks (CKM matrix) and neutrinos (PMNS
matrix) — both are determined by the same topological structure of
the internal space K = CP² × S³ × S¹.

The Hopf fibration S³ → S² is the bundle that describes the
electroweak gauge group. Its Chern-Simons invariant gives arctan(φ²)
as the natural angle associated with the golden ratio winding.

CP violation is not a free parameter in FisherGeometrics but a
topological invariant of K — it could not be otherwise.

    δ_CP = arctan(φ²) = arctan(2.618...) = 69.095°

Measured values (PDG 2024):
  CKM:  δ_CP = 69.3°  (0.30% deviation)
  PMNS: δ_CP = 69.2°  (0.14% deviation)

# Returns
CP-violating phase in degrees

# Examples
```julia
δ = cp_phase()           # → 69.095°

ckm  = ckm_wolfenstein(verbose=false)
pmns = pmns_angles(verbose=false)
ckm.δ_CP == pmns.δ_CP == cp_phase()   # → true ✓
```

See: FisherGeometrics preprint v15, sections 8 and 9.
"""
function cp_phase()
    φ = (1+sqrt(5))/2
    return rad2deg(atan(φ^2))
end
