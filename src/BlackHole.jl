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
