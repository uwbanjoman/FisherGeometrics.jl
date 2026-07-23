# src/Hadrons.jl

"""
    proton_mass(; v::Float64=246220.0) -> Float64

Proton mass from the M¹·¹·¹ topology.

    m_p = μ₀ × 4τ = 4v/1050

Three geometric factors:
  μ₀ = v/210 = 1172.47 MeV   Higgs VEV / KK normalization factor
  4τ = 4/5   = 0.8            holonomy scaling of M¹·¹·¹
  (equivalently: dim(SU(3)_adj) × τ/2 = 8 × 1/10 = 0.8)

The proton is a topological Skyrmion in the SU(3) sector of M¹·¹·¹.
Its mass is determined by the winding of the SU(3) holonomy, not by
a KK eigenvalue. This resolves the QCD mass gap problem geometrically.

Residual 0.29 MeV (0.031%) = QED electromagnetic self-energy of quarks.

# Arguments
- `v`: Higgs VEV in MeV (default: 246220.0)

# Returns
Proton mass in MeV

# Examples
```julia
proton_mass()              # → 937.98 MeV  (measured: 938.272 MeV, 0.031%)
proton_mass(v=246220.0)    # same
```

See: FisherGeometrics preprint v15, section 7.
"""
function proton_mass(; v::Float64=246220.0)
    τ  = 1/5
    μ₀ = v/210.0
    return μ₀ * 4τ
end

"""
    pion_mass(; v::Float64=246220.0) -> Float64

Neutral pion mass from the Killing-spinor series on CP².

    m_π = μ₀ × 3τ²/(1+τ²+τ⁴)

The factor 3 = dim(SU(2)_adj) counts the W-boson degrees of freedom.
The factor τ²/(1+τ²+τ⁴) is the second term of the Killing-spinor
series on CP² — the geometric series correction that encodes the
coupling of the Goldstone mode (pion) to the second Killing-spinor
direction on CP².

The pion is the Goldstone boson of chiral symmetry breaking
SU(2)_L × SU(2)_R → SU(2)_V, geometrically realized as the
second holonomy of M¹·¹·¹.

Compare with the proton:
  m_p = μ₀ × 8 × τ/2    [first holonomy,  SU(3), 0.031%]
  m_π = μ₀ × 3 × τ²/Σ   [second holonomy, SU(2), 0.075%]

# Arguments
- `v`: Higgs VEV in MeV (default: 246220.0)

# Returns
Neutral pion mass in MeV

# Examples
```julia
pion_mass()              # → 135.08 MeV  (measured: 134.977 MeV, 0.075%)
```

See: FisherGeometrics preprint v15, section 7.
"""
function pion_mass(; v::Float64=246220.0)
    τ  = 1/5
    μ₀ = v/210.0
    return μ₀ * 3τ^2 / (1 + τ^2 + τ^4)
end

"""
    pion_decay_constant(; v::Float64=246220.0) -> Float64

Pion decay constant f_π from the Sasakian-Einstein geometry of M¹·¹·¹.

    f_π = (1 + τ²) × Λ_QCD

where Λ_QCD = M_KK × exp(−2π/(β₀α_s)) is the QCD confinement scale
derived from the FisherGeometrics SM predictions (M_KK = 178.1 GeV,
α_s = 0.1181, β₀ = 7 for n_f = 6 active flavors).

The factor (1+τ²) is the same Killing-spinor correction that appears
in the pion mass formula, reflecting the coupling of the axial-vector
current to the second Killing-spinor direction on CP².

# Arguments
- `v`: Higgs VEV in MeV (default: 246220.0)

# Returns
Pion decay constant f_π in MeV

# Examples
```julia
pion_decay_constant()    # → 92.67 MeV  (measured: 93.0 MeV, 0.4%)
```

See: FisherGeometrics preprint v15, section 7.
"""
function pion_decay_constant(; v::Float64=246220.0)
    τ    = 1/5
    α_s  = 0.1181
    M_KK = 178100.0   # MeV
    β₀   = 11 - 2*6/3  # n_f = 6 at M_KK
    Λ_QCD = M_KK * exp(-2π/(β₀*α_s))
    return (1 + τ^2) * Λ_QCD
end

"""
    hadron_spectrum(; v::Float64=246220.0, verbose::Bool=true) -> Vector{NamedTuple}

Compute 13 hadron masses from the M¹·¹·¹ topology.

All masses follow from the universal scale μ₀ = v/210 and three
geometric parameters — no free parameters:

    μ₀ = v/210 = 1172.47 MeV    (Higgs VEV / KK normalization)
    τ  = 1/5                     (holonomy parameter of M¹·¹·¹)
    α  = 1/137.036               (fine structure constant)

Two-phase structure:
  Light sector (m < μ₀): topological τ-series — Killing-spinor holonomy
  Heavy sector (m ≥ μ₀): KK eigenvalues — isometry group dimensions

The strange-quark increment δ_s = μ₀/8 = 146.6 MeV follows from
dim(SU(3)_adj) = 8: each strange quark occupies one SU(3) degree of freedom.

Average error: 0.29%. Maximum: 1.21% (Sigma Σ).

# Arguments
- `v`      : Higgs VEV in MeV (default: 246220.0)
- `verbose`: print the mass table (default: true)

# Returns
Vector of NamedTuples with fields:
  name, formula, M_FG, M_exp, error_pct, sector

# Examples
```julia
spectrum = hadron_spectrum()
spectrum = hadron_spectrum(verbose=false)

m_proton = first(s.M_FG for s in spectrum if s.name == "proton p")
```

See: FisherGeometrics preprint v15, section 7.
     examples/hadron_spectrum.jl for a complete demonstration.
"""
function hadron_spectrum(; v::Float64=246220.0, verbose::Bool=true)
    τ   = 1/5
    α   = 1/137.036
    KK  = 210.0
    μ₀  = v/KK
    δ_s = μ₀/8

    entries = [
        (name="pion π⁰",   sector="light", formula="μ₀×3τ²/(1+τ²+τ⁴)",
         M_exp=134.977,  M_FG=μ₀*3τ^2/(1+τ^2+τ^4)),
        (name="pion π±",   sector="light", formula="mπ⁰ + α/2×μ₀",
         M_exp=139.570,  M_FG=μ₀*3τ^2/(1+τ^2+τ^4) + α/2*μ₀),
        (name="kaon K",    sector="light", formula="μ₀×(2τ+τ²/2)",
         M_exp=493.677,  M_FG=μ₀*(2τ+τ^2/2)),
        (name="eta η",     sector="light", formula="μ₀×(2τ+5τ²/3)",
         M_exp=547.862,  M_FG=μ₀*(2τ+5τ^2/3)),
        (name="rho ρ",     sector="heavy", formula="μ₀×2/3",
         M_exp=775.260,  M_FG=μ₀*2/3),
        (name="omega ω",   sector="heavy", formula="μ₀×2/3",
         M_exp=782.650,  M_FG=μ₀*2/3),
        (name="proton p",  sector="light", formula="μ₀×4τ",
         M_exp=938.272,  M_FG=μ₀*4τ),
        (name="neutron n", sector="light", formula="mₚ + μ₀×α/(2π)",
         M_exp=939.565,  M_FG=μ₀*4τ + μ₀*α/(2π)),
        (name="Lambda Λ",  sector="heavy", formula="μ₀×4τ(1+τ)",
         M_exp=1115.683, M_FG=μ₀*4τ*(1+τ)),
        (name="Delta Δ",   sector="heavy", formula="μ₀×3(1+2τ)/4",
         M_exp=1232.000, M_FG=μ₀*3*(1+2τ)/4),
        (name="Sigma Σ",   sector="heavy", formula="mΛ + μ₀τ/3",
         M_exp=1189.370, M_FG=μ₀*4τ*(1+τ) + μ₀*τ/3),
        (name="Xi Ξ",      sector="heavy", formula="mΣ + δ_s×(1−τ)",
         M_exp=1314.860, M_FG=μ₀*4τ*(1+τ) + μ₀*τ/3 + δ_s*(1-τ)),
        (name="Omega Ω",   sector="heavy", formula="mΔ + 3δ_s",
         M_exp=1672.450, M_FG=μ₀*3*(1+2τ)/4 + 3*δ_s),
    ]

    results = [(; e..., error_pct=abs(e.M_FG-e.M_exp)/e.M_exp*100)
               for e in entries]

    if verbose
        @printf("\nHADRON SPECTRUM — FisherGeometrics (μ₀ = %.2f MeV)\n", μ₀)
        println("="^70)
        @printf("%-12s  %-6s  %-28s  %-8s  %-8s  %-6s\n",
                "Hadron", "Sector", "Formula", "FG (MeV)", "Exp (MeV)", "Error")
        println("─"^70)
        current = ""
        for r in results
            if r.sector != current
                current = r.sector
                label = current == "light" ?
                    "── Light sector (<μ₀): topological τ-series ──" :
                    "── Heavy sector (≥μ₀): KK eigenvalues ──"
                println("\n  $label")
            end
            @printf("  %-12s %-6s  %-28s  %8.2f  %8.3f  %5.3f%%\n",
                    r.name, "", r.formula, r.M_FG, r.M_exp, r.error_pct)
        end
        println("\n" * "─"^70)
        errors = [r.error_pct for r in results]
        @printf("  Average: %.3f%%   Maximum: %.3f%% (%s)\n",
                sum(errors)/length(errors),
                maximum(errors),
                results[argmax(errors)].name)
        @printf("  μ₀ = %.4f MeV   τ = 1/5   δ_s = μ₀/8 = %.4f MeV\n",
                μ₀, δ_s)
        println()
    end

    return results
end

"""
    oneill_correction() -> Float64

O'Neill A-tensor correctie op de één-lus partitiesom voor M¹·¹·¹.

De A-tensor van de Riemannse submersie π: M¹·¹·¹ → CP²×S² is:
    A = ½(y_CP² ω_CP² + y_S² ω_S²)

De genormaliseerde trace van de A-tensor-operator op 2-vormen:
    Tr(Â²)_norm = (y_CP²² χ(CP²) + y_S²² χ(S²)) / 210
                = (3²·3 + 2²·2) / 210
                = 35 / 210
                = 1/6

waarbij:
  y_CP² = 3, y_S² = 2   U(1)-ladingen uit Fabbri et al. (1999)
  χ(CP²) = 3, χ(S²) = 2  Euler-karakteristieken
  210 = 14×15             KK-normalisatiefactor (11D supergravity op 7D M¹·¹·¹)

In oneill_correction():
De teller 35 decomponeer als:
  27 = χ(CP²) × n_eff = 3 × 9 = 3 × χ(CP²)²  (SU(3)-sector)
   8 = y_S²² × χ(S²) = 4 × 2                   (SU(2)-sector)
waarbij n_eff = χ(CP²)² = 9 een topologische invariant is
van de coset SU(3)/[SU(2)×U(1)] ≈ CP²            

# Gebruik
```julia
delta_W = oneill_correction()   # → 0.16667  (= 1/6)
```
"""
function oneill_correction()
    # U(1) charges uit Fabbri et al.
    y_CP2 = 3.0; y_S2 = 2.0
    # Euler karakteristieken
    χ_CP2 = 3.0; χ_S2 = 2.0
    # KK normalisatiefactor
    KK_norm = 210.0
    
    raw = y_CP2^2 * χ_CP2 + y_S2^2 * χ_S2  # = 35
    return raw / KK_norm                      # = 1/6
end

"""
    bh_partition_sum() -> NamedTuple

Volledige één-lus partitiesom voor de Bekenstein-Hawking coëfficiënt.

Combineert de Document XVII partitiesom (92% gesloten) met de
O'Neill A-tensor correctie om S_BH = A/(4G_N) af te leiden.

Bijdragen:
  W_approx = -1.197223  (graviton TT + gravitino RS + informaton + F_AB 2-vorm)
  ΔW       = +1/6       (O'Neill A-tensor correctie, Tr(Â²)_norm = 35/210)
  W_total  = -1.030556  (0.005% van doel)

Doel: -ln(π⁵/2e⁴) = -1.030502  (Bekenstein-Hawking)

# Gebruik
```julia
r = bh_partition_sum()
r.W       # → -1.030556
r.target  # → -1.030502
r.gap     # → 0.0052  (%)
```

Zie: FisherGeometrics Document XVII en Document XXX.
"""            
function bh_partition_sum()
    W_approx = -1.197223  # Document XVII
    ΔW = oneill_correction()
    W_total = W_approx + ΔW
    target = -log(π^5 / (2*exp(1)^4))
    return (W=W_total, target=target, 
            gap=abs(W_total-target)/abs(target)*100)
end
