# src/States.jl

"""
    rho_KK(R; n=6) -> Matrix{ComplexF64}

KK-thermische dichtheidsmatrix bij compactificatieschaal R:
    ρ(R) = diag(1, e^{-16/R}, e^{-16/R}, e^{-16/R}, e^{-144/R}, e^{-144/R}) / Z

R→0: ρ → |0⟩⟨0|  (massaloos gravitino)
R→∞: ρ → I/6     (vacuüm ρ*)
"""
function rho_KK(R::Real; n::Int=6)
    p0=1.0; p1=exp(-16/R); p2=exp(-144/R); Z=p0+3*p1+2*p2
    return Matrix{ComplexF64}(Diagonal([p0,p1,p1,p1,p2,p2]./Z))
end

"""
    rho_KK_eigenvalues(R; n=6) -> Vector{Float64}

Geeft de eigenwaarden van ρ(R) als vector terug (zonder matrix-constructie).
Handig voor analytische berekeningen.
"""
function rho_KK_eigenvalues(R::Real; n::Int=6)
    p0=1.0; p1=exp(-16/R); p2=exp(-144/R); Z=p0+3*p1+2*p2
    return [p0,p1,p1,p1,p2,p2]./Z
end

"""
    D2_bures_KK(R; n=6) -> Float64

Exacte analytische Bures-afstand tussen ρ(R) en ρ* = I/6:

    D²(R) = 2(1 − (1+3e^{-8/R}+2e^{-72/R}) / √(6Z))

Geldig voor alle R > 0.
"""
function D2_bures_KK(R::Real; n::Int=6)
    p1=exp(-16/R); p2=exp(-144/R); Z=1+3*p1+2*p2
    return max(2*(1-(1/sqrt(n*Z))*(1+3*exp(-8/R)+2*exp(-72/R))), 0.0)
end

"""
    rho_complex_KK(R, phi) -> Matrix{ComplexF64}

KK-toestand met een complexe fase-draai in de (1,2)-component.
Wordt gebruikt voor het testen van off-diagonale perturbaties.
"""
function rho_complex_KK(R::Real, phi::Real)
    p0=1.0; p1=exp(-16/R); p2=exp(-144/R); Z=p0+3*p1+2*p2
    ρ = Matrix{ComplexF64}(Diagonal([p0,p1,p1,p1,p2,p2]./Z))
    val = 0.05 * exp(im * phi)
    ρ[1,2] = val; ρ[2,1] = conj(val)
    return ρ / tr(ρ)
end

"""
    entropy_KK(R; n=6) -> Float64

Von Neumann entropie van de KK-familie ρ(R), direct uit de eigenwaarden.
Sneller dan entropy(rho_KK(R)) voor grote berekeningen.
"""
function entropy_KK(R::Real; n::Int=6)
    p = rho_KK_eigenvalues(R; n=n)
    return -sum(v * log(v) for v in p if v > 1e-15)
end

"""
    kk_path(R_start, R_end; n_steps=10) -> Vector{Matrix{ComplexF64}}

Geeft een reeks KK-toestanden terug langs het pad R_start → R_end.
Handig voor het visualiseren van de BGK-relaxatie.
"""
function kk_path(R_start::Real, R_end::Real; n_steps::Int=10)
    Rs = range(R_start, R_end, length=n_steps)
    return [rho_KK(R) for R in Rs]
end
