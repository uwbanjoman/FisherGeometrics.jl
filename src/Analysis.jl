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
