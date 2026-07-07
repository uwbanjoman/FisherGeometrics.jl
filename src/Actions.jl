# src/Actions.jl

function information_action(g::FisherMetric,
                            rhos,
                            basis,
                            L;
                            Δ = 1.0)

    I = 0.0

    for ρ in rhos

        G   = metric_matrix(g,ρ,basis)
        Ric = ricci(g,ρ,basis)
        S   = scalar_curvature(g,ρ,basis)

        I += L(ρ,G,Ric,S) * Δ

    end

    return I

end
