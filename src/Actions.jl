# src/Actions.jl

function information_action(g::FisherMetric,
                            rhos,
                            basis,
                            L;
                            Δ = 1.0)
    I = 0.0

    for ρ in rhos
        # 1. Bereken de metriekmatrix
        G = metric_matrix(g, ρ, basis)
        
        # 2. Volume-element (invariantie onder coördinatentransformaties)
        # Gebruik abs() voor de Riemannse (positief-definiete) Fisher-metriek
        det_G = det(G)
        dV = sqrt(abs(det_G)) * Δ

        # 3. Bereken de krommingen 
        # Tip: Als je pakket dit ondersteunt, bereken Ric en S in één pass 
        Ric = ricci(g, ρ, basis)
        S   = scalar_curvature(g, ρ, basis)

        # 4. Voeg toe aan de totale actie, inclusief de invariante maat
        I += L(ρ, G, Ric, S) * dV
    end

    return I
end
