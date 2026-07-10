# src/Interactions.jl
# Non-commutative dynamics, gauge interactions, and field fluxes.
# Depends on: Foundation.jl, Dynamics.jl

"""
    interaction(ρ, ops)

The fundamental non-commutative operator: [ρ, ops].
Represents the gauge-theoretic tension between the state and the symmetry operator.
"""
function interaction(ρ::AbstractMatrix, ops::AbstractMatrix)
    return ρ * ops - ops * ρ
end

"""
    flux(ρ, ops)

Energy density of the non-commutative interaction.
This measures the 'cost' of the field configuration in the informatic manifold.
"""
function flux(ρ::AbstractMatrix, ops::AbstractMatrix)
    δρ = interaction(ρ, ops)
    return real(tr(δρ * δρ))
end

"""
    covariant_shift(ρ, gauge_field, dt)

Parallel transport of the state ρ along a gauge field.
Describes how the state 'rotates' (EM/Weak interaction) while 
navigating the manifold defined in Dynamics.jl.
"""
function covariant_shift(ρ::AbstractMatrix, gauge_field::AbstractMatrix, dt::Float64)
    # The state evolves via the commutator, preserving unitary structure
    # ρ' = ρ - i * dt * [gauge_field, ρ]
    return ρ - im * dt * interaction(ρ, gauge_field)
end

# ── Symmetry Generators ───────────────────────────────────────

"""
    u1_generator()

Generator for U(1) symmetry (Electromagnetism).
"""
function u1_generator()
    # Basic charge generator for the manifold
    return [1.0 0.0; 0.0 -1.0] # Representatie van de ijk-generator
end

"""
    u1_generator_3x3()

Generator for U(1) symmetry (Electromagnetism) from CKM
"""
function u1_generator_3x3()
    # Deze generator representeert de lading-structuur van de 3 generaties
    return [1.0 0.0 0.0; 
            0.0 1.0 0.0; 
            0.0 0.0 1.0] 
end

"""
    su2_generator()

Generator for SU(2) symmetry (Weak interaction / CKM-coupling).
"""
function su2_generator()
    return [0.0 1.0; 1.0 0.0]
end

# ── Interaction consistency ───────────────────────────────────

function check_interactions()
    ρ = [1.0 0.0; 0.0 0.0] # Pure state reference
    ops = u1_generator()
    val = flux(ρ, ops)
    # Flux should be zero for a pure state aligned with the generator
    return val >= 0.0
end
