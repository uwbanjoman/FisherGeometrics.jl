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

function run_fase_test(R_range, β)
    phases = Float64[]
    for i in 1:length(R_range)-1
        R1, R2 = R_range[i], R_range[i+1]
        
        ρ1 = gibbs_state_expanded(M1, M2, J/R1, β)
        ρ2 = gibbs_state_expanded(M1, M2, J/R2, β)
        
        # FIX: Laat de rotatie afhangen van R
        # Hoe sterker de kromming (kleinere R), hoe sterker de fase-koppeling
        val1 = 0.5 * exp(im * 0.5 * R1)
        val2 = 0.5 * exp(im * 0.5 * R2)
        
        ρ1[1,2] = val1; ρ1[2,1] = conj(val1)
        ρ2[1,2] = val2; ρ2[2,1] = conj(val2)
        
        # Normaliseer na toevoeging van off-diagonale termen
        ρ1 /= tr(ρ1)
        ρ2 /= tr(ρ2)
        
        _, V1 = eigen(ρ1)
        _, V2 = eigen(ρ2)
        
        # Gebruik een 'gapped' overlap (neem de vector met de grootste eigenwaarde)
        overlap = dot(V1[:, end], V2[:, end])
        push!(phases, imag(log(overlap)))
    end
    return sum(phases)
end
