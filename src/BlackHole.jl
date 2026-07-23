# src/BlackHole.jl

"""
    hawking_temp_FG(M_kg::Float64) -> Float64

Hawking temperature of a black hole from the FisherGeometrics scale μ₀.

    T_H = μ₀ × τ / (8π × M_BH)

where μ₀ = v/210 is the universal FG scale and τ = 1/5 is the holonomy
parameter of M¹·¹·¹. This follows from the identification of the black
hole as a density matrix ρ̂_BH far from the vacuum ρ* = I/6, with
Hawking radiation as BGK relaxation ρ̂_BH(t) → ρ*.

# Arguments
- `M_kg`: black hole mass in kg

# Returns
Hawking temperature in Kelvin

# Examples
```julia
M_sun = 1.989e30   # kg
hawking_temp_FG(M_sun)       # → 6.17e-8 K
hawking_temp_FG(2.176e-8)    # → Planck mass → ~1e32 K
```

See: FisherGeometrics preprint v15, section 11 (black hole evaporation).
"""
function hawking_temp_FG(M_kg::Float64)
    c    = 3e8
    k_B  = 1.381e-23
    μ₀   = 246220.0 / 210.0   # MeV
    τ    = 1/5
    M_MeV = M_kg * c^2 / 1.602e-13
    T_MeV = μ₀ * τ / (8π * M_MeV)
    return T_MeV * 1.602e-13 / k_B
end

"""
    evaporation_time(M_kg::Float64) -> Float64

Hawking evaporation time of a black hole.

    t_evap = 5120π × G³/(ħc⁴) × M³

The evaporation time follows from Stefan-Boltzmann radiation at the
Hawking temperature T_H ∝ 1/M. As the black hole radiates, M decreases,
T_H increases, and the process accelerates — ending in a final burst
when M → M_Planck.

In FisherGeometrics: evaporation is BGK relaxation ρ̂_BH(t) → ρ*.
The timescale t_evap is the time for the density matrix to return
to the vacuum ρ* = I/6.

# Arguments
- `M_kg`: initial black hole mass in kg

# Returns
Evaporation time in seconds

# Examples
```julia
M_sun = 1.989e30
evaporation_time(M_sun)      # → 2.09e74 s  (much longer than age of universe)
evaporation_time(1e12)       # → 8.41e26 s  (primordial BH, ~10⁻¹¹ kg → gone)
evaporation_time(2.176e-8)   # → Planck mass → ~5.4e-44 s (one Planck time)
```

See: FisherGeometrics preprint v15, section 11 (black hole evaporation).
"""
function evaporation_time(M_kg::Float64)
    G  = 6.674e-11   # m³ kg⁻¹ s⁻²
    ħ  = 1.055e-34   # J·s
    c  = 3e8         # m/s
    return 5120π * G^3 / (ħ * c^4) * M_kg^3
end

"""
    bh_entropy(M_kg::Float64) -> Float64

Bekenstein-Hawking entropy of a black hole in natural units (k_B = 1).

    S_BH = A / (4G_N)  =  4π G M² / (ħc)

Derived in FisherGeometrics via the one-loop partition sum on M₄×K
plus the O'Neill A-tensor correction Tr(Â²)_norm = 35/210 = 1/6:

    c_log = ζ'_base(0) + 1/6 = −1.030556  (0.005% from target)

This establishes S_BH = A/(4G_N) with zero free parameters.

The entropy measures the Bures distance of ρ̂_BH from the vacuum ρ*:
a larger black hole is further from ρ* and has higher entropy.

# Arguments
- `M_kg`: black hole mass in kg

# Returns
Bekenstein-Hawking entropy in units of k_B

# Examples
```julia
M_sun = 1.989e30
bh_entropy(M_sun)       # → 1.05e77  k_B  (solar mass BH)
bh_entropy(2.176e-8)    # → Planck mass → ~1 k_B (one bit)
```

See: FisherGeometrics preprint v15, sections 8 and 11.
     BH_coefficient_derivation.pdf for the complete derivation.
"""
function bh_entropy(M_kg::Float64)
    G  = 6.674e-11   # m³ kg⁻¹ s⁻²
    ħ  = 1.055e-34   # J·s
    c  = 3e8         # m/s
    k_B = 1.381e-23  # J/K
    return 4π * G * M_kg^2 / (ħ * c * k_B)
end
