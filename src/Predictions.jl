# src/Predictions.jl

"""
    glueball_mass(; v::Float64=246220.0) -> Float64

Predicted mass of the scalar glueball 0⁺⁺ from M¹·¹·¹ geometry.

    m_G = μ₀ × √2

The scalar glueball is a pure eigenmode of the CP² metric without
quark winding. The √2 factor arises from the S²/CP² radius ratio
in the Freund-Rubin compactification (R²_S² = R²_CP²/4), the same
geometric ratio that gives the Regge slope α' = (2/3)μ₀².

The f₀(1710) resonance at 1710 MeV is the primary experimental
candidate. The 52 MeV difference (3.0%) is attributable to mixing
with nearby quark-antiquark states. A pure glueball state should
be sought at 1658 MeV in gluon-rich decay processes (J/ψ → γ + X,
pp̄ → ggg → hadrons).

# Arguments
- `v`: Higgs VEV in MeV (default: 246220.0)

# Returns
Predicted glueball mass in MeV

# Examples
```julia
glueball_mass()     # → 1658.13 MeV  (f₀(1710): 1710 MeV, 3.0%)
```

See: FisherGeometrics preprint v15, section 10 (experimental predictions).
"""
function glueball_mass(; v::Float64=246220.0)
    μ₀ = v/210.0
    return μ₀ * sqrt(2)
end

"""
    dark_baryon_mass(; v::Float64=246220.0) -> Float64

Predicted mass of the dark baryon — a falsifiable FisherGeometrics prediction.

    m_DM = μ₀ × τ = v/(210 × 5) = v/1050

The dark baryon is a topologically stable neutral particle stabilized
by the SU(3) topological index of M¹·¹·¹. Unlike ordinary baryons
(which wind around the SU(3) sector with factor 4τ), the dark baryon
occupies the fundamental τ winding — the minimal stable configuration.

Properties:
  Mass:     234.5 MeV  (between pion and kaon)
  Charge:   0          (no U(1) coupling — invisible to photons)
  Stability: SU(3) topological index (cannot decay to lighter hadrons)
  Coupling: gravitational only + rare strong-sector interactions

This makes it a cold dark matter candidate detectable via:
  - Gravitational effects (galaxy rotation curves)
  - Rare hadronic decays in high-energy colliders
  - Direct detection via nuclear recoil (mass ~ 234.5 MeV)

Current experimental status: unknown — falsifiable prediction.
Search channels: BaBar, Belle II, LHCb missing-mass searches.

# Arguments
- `v`: Higgs VEV in MeV (default: 246220.0)

# Returns
Predicted dark baryon mass in MeV

# Examples
```julia
dark_baryon_mass()     # → 234.50 MeV  (unknown experimentally)
```

See: FisherGeometrics preprint v15, section 10 (experimental predictions).
"""
function dark_baryon_mass(; v::Float64=246220.0)
    τ  = 1/5
    μ₀ = v/210.0
    return μ₀ * τ
end

"""
    regge_slope(; v::Float64=246220.0) -> Float64

Universal Regge slope α' from the S²/CP² radius ratio in M¹·¹·¹.

    α' = (2/3) × μ₀²  =  4 × (μ₀²/6)

Two geometric factors:
  1/6  =  Tr(Â²)_norm = 35/210  O'Neill A-tensor torsion of M¹·¹·¹
  ×4   =  R²_CP²/R²_S²          S²/CP² radius ratio (Freund-Rubin)

The factor 4 = R²_CP²/R²_S² follows from the Freund-Rubin
compactification of 11D supergravity on M¹·¹·¹: the S² fiber has
radius R_S² = R_CP²/2, so mass-squared excitations on S² are enhanced
by a factor 4 relative to the CP² base scale.

Together: α' = 4 × μ₀²/6 = (2/3)μ₀²

This connects the string tension of the QCD flux tube directly
to the 7D geometry of M¹·¹·¹ — no free parameters.

Hadronic Regge trajectories satisfy:
    m²_n = m²_0 + n × α'

# Arguments
- `v`: Higgs VEV in MeV (default: 246220.0)

# Returns
Regge slope α' in GeV²

# Examples
```julia
regge_slope()     # → 0.9165 GeV²  (experimental: 0.9–1.1 GeV², 1.8%)
```

See: FisherGeometrics preprint v15, sections 7.4 and 10.
"""
function regge_slope(; v::Float64=246220.0)
    μ₀ = v/210.0   # MeV
    return (2/3) * (μ₀/1000)^2   # GeV²
end

"""
    mass_gap_floor(; v::Float64=246220.0) -> Float64

Fundamental energy threshold for non-Goldstone strong resonances.

    E_gap = μ₀/6

No strong resonance (except Goldstone bosons such as the pion, kaon,
and eta) can exist below this threshold. The factor 1/6 = Tr(Â²)_norm
is the O'Neill A-tensor torsion of M¹·¹·¹ — the same factor that
closes the Bekenstein-Hawking gap.

This explains the observed gap in the hadron spectrum between the
pion (135 MeV) and the kaon (494 MeV): no non-Goldstone resonance
exists below E_gap = 195 MeV.

    E_gap = μ₀/6 = 195.4 MeV

The pion (135 MeV) lies below E_gap because it is a Goldstone boson
of chiral symmetry breaking — protected by topology, not by the
KK floor. All other hadrons lie above E_gap.

Consistency check:
  pion  π⁰: 135 MeV < E_gap  (Goldstone boson — exempt) ✓
  kaon   K: 494 MeV > E_gap  ✓
  rho    ρ: 775 MeV > E_gap  ✓
  proton p: 938 MeV > E_gap  ✓

# Arguments
- `v`: Higgs VEV in MeV (default: 246220.0)

# Returns
Mass gap floor E_gap in MeV

# Examples
```julia
mass_gap_floor()    # → 195.41 MeV

# Verify: all non-Goldstone hadrons lie above E_gap
E_gap = mass_gap_floor()
@printf("Kaon (494 MeV) > E_gap (%.1f MeV): %s\\n",
        E_gap, 493.677 > E_gap ? "✓" : "✗")
```

See: FisherGeometrics preprint v15, section 10 (experimental predictions).
"""
function mass_gap_floor(; v::Float64=246220.0)
    μ₀ = v/210.0
    return μ₀/6
end
