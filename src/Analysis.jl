# src/Analysis.jl

"""
    compute_I_R_exact(R_values; n=6)

Berekent de Fisher-informatie ratio I(R) voor een gegeven bereik van R.
Geeft een vector van NamedTuples terug voor verdere analyse/plotting.
"""
function compute_I_R_exact(R_values::Vector{<:Real}; n::Int=6)
    T = su_basis(n)
    results = Vector{NamedTuple}(undef, length(R_values))
    
    for (i, R) in enumerate(R_values)
        ρ  = rho_KK_exact(R; n=n)
        D2 = D2_KK(R; n=n)
        SF = SF_exact(ρ, T)
        
        # Voorkom deling door nul en definieer I
        I_R = D2 > 1e-12 ? SF / D2 : Inf
        
        results[i] = (R=R, D2=D2, SF=SF, I=I_R)
    end
    
    return results
end

"""
    print_I_R_summary(results)

Verantwoordelijk voor de console-output. Zo blijft je berekeningsfunctie puur.
"""
function print_I_R_summary(results::Vector{<:NamedTuple})
    println("═"^60)
    println("Fisher Informatie Analyse: S_F(ρ) / D²")
    println("═"^60)
    
    for r in results
        @printf("R = %8.1f | SF = %8.4f | D² = %8.5f | I = %8.3f\n", 
                r.R, r.SF, r.D2, r.I)
    end
    
    I_vals = [r.I for r in results]
    idx = argmin(I_vals)
    println("─"^60)
    @printf("Minimum I(R) bij R* ≈ %.1f (I = %.3f)\n", results[idx].R, results[idx].I)
end

"""
    run_fase_test(R_range, β, M1, M2, J)

Calculates the total geometric (Berry) phase of a Kaluza-Klein black hole model 
as the system traverses a radial path defined by `R_range`.

This function simulates the interaction between the thermal Gibbs state of the 
su(2) subsystem and the geometric curvature of spacetime. The phase is determined 
by summing the local phase shifts (holonomy) of the principal eigenvector of the 
density matrix ρ(R).

# Arguments
- `R_range`: A vector of radial coordinates along which the phase is integrated.
- `β`: The inverse temperature (1/T).
- `M1`, `M2`: Kaluza-Klein mass modes defining the energy scale.
- `J`: The coupling constant or angular momentum of the system.

# Returns
- `total_phase`: The accumulated geometric phase in radians. A non-zero value 
  indicates the presence of topological charge or curvature within the event horizon.

# Example
```julia
phases = run_fase_test(exp.(range(log(50), log(4260), length=50)), 2.0, 1, 1, 0.5)
"""
function run_fase_test(R_range, β, M1, M2, J)
    total_phase = 0.0
    for i in 1:length(R_range)-1
        R1, R2 = R_range[i], R_range[i+1]
        
        # Bereken je toestanden
        ρ1 = gibbs_state_expanded(M1, M2, J/R1, β)
        ρ2 = gibbs_state_expanded(M1, M2, J/R2, β)
        
        # Voeg de R-afhankelijke fase toe
        val1 = 0.5 * exp(im * (2π * R1 / 4260)) 
        val2 = 0.5 * exp(im * (2π * R2 / 4260))
        ρ1[1,2] += val1; ρ1[2,1] += conj(val1)
        ρ2[1,2] += val2; ρ2[2,1] += conj(val2)
        ρ1 /= tr(ρ1); ρ2 /= tr(ρ2)
        
        # Eigenvectoren
        _, V1 = eigen(ρ1)
        _, V2 = eigen(ρ2)
        
        # Overlap
        overlap = dot(V1[:, end], V2[:, end])
        
        # De 'Berry-fase' stap is het verschil in hoek
        # Dit vangt de verandering op in plaats van de absolute waarde
        step_phase = angle(overlap)
        total_phase += step_phase
    end
    return total_phase
end
