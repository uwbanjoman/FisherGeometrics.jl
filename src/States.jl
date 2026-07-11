# src/States.jl

"""
    rho_KK(R; n=6) -> Matrix{ComplexF64}

KK-thermische dichtheidsmatrix bij compactificatieschaal R:
    ρ(R) = diag(1, e^{-16/R}, ..., e^{-144/R}) / Z
R→0: ρ → |0><0|  (massaloos gravitino)
R→∞: ρ → I/6     (vacuüm ρ*)
"""
function rho_KK(R::Real; n::Int=6)
    p0=1.0; p1=exp(-16/R); p2=exp(-144/R); Z=p0+3*p1+2*p2
    return Matrix{ComplexF64}(Diagonal([p0,p1,p1,p1,p2,p2]./Z))
end

function entropy(ρ::Diagonal)
    p = ρ.diag
    return -sum(x > 0 ? x * log(x) : 0.0 for x in p)
end
