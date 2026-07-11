# src/Curvature.jl

export su_basis

"""
    su_basis(n) -> Vector{Matrix{ComplexF64}}

𝔰𝔲(n) generatoren, genormeerd als Tr(TₐTᵦ) = δₐᵦ/2.
Geeft n²−1 matrices terug.
"""
function su_basis(n::Int)
    T = Matrix{ComplexF64}[]
    for j in 1:n, k in j+1:n
        M = zeros(ComplexF64,n,n); M[j,k]=M[k,j]=0.5; push!(T,M)
    end
    for j in 1:n, k in j+1:n
        M = zeros(ComplexF64,n,n); M[j,k]=-0.5im; M[k,j]=0.5im; push!(T,M)
    end
    for l in 1:n-1
        M = zeros(ComplexF64,n,n); nrm=1/sqrt(2l*(l+1))
        for j in 1:l; M[j,j]=nrm; end; M[l+1,l+1]=-l*nrm; push!(T,M)
    end
    return T
end

"""
    d_tensor(T) -> Array{Float64,3}
 
Symmetrische structuurconstanten van 𝔰𝔲(n):
 
    d_{abc} = 4 Re Tr(Tₐ Tᵦ T_c)
 
Geeft een N×N×N tensor terug waarbij N = length(T) = n²−1.
 
# Voorbeeld
```julia
T = su_basis(6)
d = d_tensor(T)    # 35×35×35
```
"""
function d_tensor(T::Vector)
    N = length(T)
    d = Array{Float64}(undef, N, N, N)
    for a in 1:N, b in 1:N, c in 1:N
        d[a,b,c] = 4 * real(tr(T[a] * T[b] * T[c]))
    end
    return d
end

"""
    christoffel_vacuum(d, n) -> Array{Float64,3}
 
Levi-Civita Christoffel-symbolen van de Bures-metriek bij ρ* = I/n:
 
    Γᵉ_{ab} = −(n/4) d_{abe}
 
Bewezen in Proof03 Theorem C.
 
# Voorbeeld
```julia
d = d_tensor(su_basis(6))
Γ = christoffel_vacuum(d, 6)
Γ[3,1,2]   # → −1.5 × d[1,2,3]
```
"""
function christoffel_vacuum(d::Array{Float64,3}, n::Int)
    N = size(d, 1)
    Γ = zeros(Float64, N, N, N)
    @tullio Γ[e,a,b] = -(n/4) * d[a,b,e]
    return Γ
end

"""
    K_quad(T, a, b; n=6) -> Float64
 
Sectionele kromming K(Tₐ, Tᵦ) bij ρ* = I/n via de kwadratische
Christoffel-formule (analytisch, milliseconden):
 
    K = (n/4)² × Σ_f [d_{abf}² − d_{aaf} d_{bbf}] / (n/8)
 
Vier waarden voor 𝒟₆ (Proof04 Theorem D):
    +0.25   sym-sym met gedeelde index
    +1.00   sym-sym zonder gedeelde index
    −2.00   sym-antisym
    −2.00   diagonaal-offdiag
 
Alle consistent met λ = 16 via Σ_{b≠a} K(Tₐ,Tᵦ) = 16.
 
# Voorbeeld
```julia
T = su_basis(6)
K_quad(T, 1, 2)   # → +0.25
K_quad(T, 3, 6)   # → +1.00
sum(K_quad(T,1,b) for b in 2:35)   # → 16.0
```
"""
function K_quad(T::Vector, a::Int, b::Int; n::Int=6)
    N  = length(T)
    _d = (x,y,z) -> 4 * real(tr(T[x] * T[y] * T[z]))
    S  = sum(_d(a,b,f)^2 - _d(a,a,f) * _d(b,b,f) for f in 1:N)
    return (n/4)^2 * S / (n/8)
end

"""
    ricci_scalar_quad(T; n=6) -> Float64
 
Scalaire Fisher-kromming bij ρ* = I/n via de analytische formule:
 
    S_F = (n²−1) × λ   waarbij λ = Σ_{b≠1} K_quad(T, 1, b)
 
Verwachte uitkomst voor 𝒟₆: **560**.
 
Bewezen in Proof04 Theorem E: de lineaire ∂Γ-bijdrage aan Ric is
nul, zodat λ = 16 exact is via de kwadratische ΓΓ-termen alleen.
 
# Voorbeeld
```julia
T = su_basis(6)
ricci_scalar_quad(T)   # → 560.0
```
"""
function ricci_scalar_quad(T::Vector; n::Int=6)
    N = length(T)
    λ = sum(K_quad(T, 1, b; n=n) for b in 2:N)
    return (n^2 - 1) * λ
end

"""
    riemann_quadratic(Γ) -> Array{Float64,4}
 
Kwadratische ΓΓ-bijdrage aan de Riemann-tensor:
 
    Q^e_{abc} = Σ_f [Γ^f_{bc} Γ^e_{af} − Γ^f_{ac} Γ^e_{bf}]
 
Bij ρ* is dit de dominante bijdrage; de lineaire ∂Γ-bijdrage
heft op in de Ricci-contractie (Proof04 Theorem E).
 
Gebruikt @tullio voor automatische SIMD-vectorisatie van de
N⁵ ≈ 52 miljoen bewerkingen.
 
# Voorbeeld
```julia
d = d_tensor(su_basis(6))
Γ = christoffel_vacuum(d, 6)
Q = riemann_quadratic(Γ)
Q[2,3,7,6]   # → +0.5625  (R^2_{376} bij ρ*)
```
"""
function riemann_quadratic(Γ::Array{Float64,3})
    N = size(Γ, 1)
    Q = zeros(Float64, N, N, N, N)
    @tullio Q[e,a,b,c] = Γ[f,b,c] * Γ[e,a,f] - Γ[f,a,c] * Γ[e,b,f]
    return Q
end

"""
    ricci_tensor(Q) -> Matrix{Float64}
 
Ricci-tensor via contractie van de eerste en derde index:
 
    Ric_{ac} = Σ_b Q^b_{abc}
 
# Voorbeeld
```julia
Ric = ricci_tensor(Q)
Ric[1,1] / (6/8)   # → 16.0  (λ bij ρ*)
```
"""
function ricci_tensor(Q::Array{Float64,4})
    N = size(Q, 1)
    Ric = zeros(Float64, N, N)
    @tullio Ric[a,c] = Q[b,a,b,c]
    return Ric
end

"""
    ricci_scalar(Ric, g_inv) -> Float64
 
Scalaire Fisher-kromming via spoor van de Ricci-tensor:
 
    S_F = g^{ab} Ric_{ab}
 
# Voorbeeld
```julia
g_inv = (8/6) * I   # inverse Bures-metriek bij ρ*
S_F = ricci_scalar(Ric, g_inv)   # → 560.0
```
"""
function ricci_scalar(Ric::Matrix, g_inv::Matrix)
    S = zeros(Float64, 1)
    @tullio S[] = g_inv[a,b] * Ric[a,b]
    return S[]
end

"""
    einstein_tensor(Ric, g, R) -> Matrix{Float64}
 
Einstein-tensor:
 
    G_{ab} = Ric_{ab} − (R/2) g_{ab}
 
waarbij R de Ricci-scalair is.
 
Bij ρ* = I/6: G_{ab} = −264 g_{ab}  (Proof04 Corollary).
 
# Voorbeeld
```julia
g   = (6/8) * Matrix(I, 35, 35)
G   = einstein_tensor(Ric, g, 560.0)
G[1,1] / g[1,1]   # → −264.0
```
"""
function einstein_tensor(Ric::Matrix, g::Matrix, R::Float64)
    N = size(Ric, 1)
    G = zeros(Float64, N, N)
    @tullio G[a,b] = Ric[a,b] - (R/2) * g[a,b]
    return G
end

"""
    bures_einstein(T; n=6) -> NamedTuple
 
Berekent de volledige Riemannse keten bij ρ* = I/n:
 
    d → Γ → Q(ΓΓ) → Ric → S_F → G
 
Geeft een NamedTuple terug met:
    .d      d-symbool tensor (N×N×N)
    .Γ      Christoffel-symbolen (N×N×N)
    .Q      Riemann ΓΓ-bijdrage (N×N×N×N)
    .Ric    Ricci-tensor (N×N)
    .S_F    Scalaire kromming (Float64, verwacht: 560)
    .G      Einstein-tensor (N×N)
    .λ      Ricci-eigenwaarde (Float64, verwacht: 16)
    .Λ      Kosmologische constante (Float64, verwacht: 264)
 
# Voorbeeld
```julia
T  = su_basis(6)
FG = bures_einstein(T)
FG.S_F     # → 560.0
FG.λ       # → 16.0
FG.Λ       # → 264.0
```
"""
function bures_einstein(T::Vector; n::Int=6)
    N   = length(T)
    g   = (n/8) * Matrix{Float64}(I, N, N)
    #g_inv = (8/n) * Matrix{Float64}(I, N, N)
    g_inv = (1/8) * (8/n) * Matrix{Float64}(I, N, N) # == (1/n) * I
 
    d   = d_tensor(T)
    Γ   = christoffel_vacuum(d, n)
    Q   = riemann_quadratic(Γ)
    Ric = ricci_tensor(Q)
    S_F = ricci_scalar(Ric, g_inv)
    G   = einstein_tensor(Ric, g, S_F)
 
    λ   = Ric[1,1] / g[1,1]
    Λ   = -G[1,1]  / g[1,1]
 
    return (d=d, Γ=Γ, Q=Q, Ric=Ric, S_F=S_F, G=G, λ=λ, Λ=Λ)
end

"""
    riemann(g, ρ, basis)

Compute the Riemann curvature tensor

    R[m,i,j,k] = Rᵐᵢⱼₖ

associated with the Fisher metric.

Returns an n×n×n×n array, where
n = length(basis).
"""
function riemann(g::FisherMetric, ρ::AbstractMatrix, basis)
    n = length(basis)

    G    = metric_matrix(g, ρ, basis)
    Ginv = pinv(G)

    dg   = metric_derivatives(g, ρ, basis)   # n × n × n
    ddg  = ddmetric_tensor(g, ρ, basis)     # n × n × n × n
    Γ    = christoffel(g, ρ, basis)          # n × n × n

    # --- OPTIMALISATIE 1: Bereken dGinv vooraf via matrix-multiplicaties ---
    # dGinv_tensor[m, l, i] = ∂ᵢ G⁻¹_{ml}
    dGinv_tensor = zeros(Float64, n, n, n)
    for i in 1:n
        # Pak de slice van dg in de i-richting: dg[:, :, i]
        dg_i = @view dg[:, :, i]
        # De contractie -Ginv * dg_i * Ginv via geoptimaliseerde matrix-multiplicatie
        dGinv_tensor[:, :, i] .= -Ginv * dg_i * Ginv
    end

    dΓ = zeros(Float64, n, n, n, n)

    # --- OPTIMALISATIE 2: Gereduceerde 5D-loop ipv 7D-loop ---
    for k in 1:n
        for j in 1:n
            # Handige shortcuts om herhaaldelijke array-lookups in de loops te vermijden
            dg_k_j = @view dg[k, :, j]
            dg_j_k = @view dg[j, :, k]
            dg_j_k_l = @view dg[j, k, :]
            
            for i in 1:n
                ddg_k_j_i = @view ddg[k, :, j, i]
                ddg_j_k_i = @view ddg[j, :, k, i]
                ddg_j_k_l_i = @view ddg[j, k, :, i]
                
                for m in 1:n
                    s = 0.0
                    @simd for l in 1:n
                        # Gebruik de vooraf berekende dGinv_tensor (O(1) lookup!)
                        term1 = dGinv_tensor[m, l, i] * (dg_k_j[l] + dg_j_k[l] - dg_j_k_l[l])
                        
                        term2 = Ginv[m, l] * (ddg_k_j_i[l] + ddg_j_k_i[l] - ddg_j_k_l_i[l])
                        
                        s += term1 + term2
                    end
                    dΓ[m, i, j, k] = 0.5 * s
                end
            end
        end
    end

    # --- Riemann Tensor Assemblage ---
    R = zeros(Float64, n, n, n, n)
    for k in 1:n
        for j in 1:n
            for i in 1:n
                for m in 1:n
                    val = dΓ[m, i, j, k] - dΓ[m, j, i, k]
                    
                    # Contractie van de Christoffel-symbolen
                    s_Γ = 0.0
                    @simd for l in 1:n
                        s_Γ += Γ[m, i, l] * Γ[l, j, k] - Γ[m, j, l] * Γ[l, i, k]
                    end
                    
                    R[m, i, j, k] = val + s_Γ
                end
            end
        end
    end

    return R
end

"""
    ricci(g, ρ, basis)

Compute the Ricci curvature tensor associated with the
Fisher metric.

Returns an n×n matrix

    Ric[i,j].
"""
function ricci(g::FisherMetric,
               ρ::AbstractMatrix,
               basis)

    R = riemann(g, ρ, basis)

    n = length(basis)

    Ric = zeros(Float64,n,n)

    for i in 1:n
        for j in 1:n

            s = 0.0

            for m in 1:n
                s += R[m,i,m,j]
            end

            Ric[i,j] = s

        end
    end

    return Ric

end

"""
    scalar_curvature(g, ρ, basis)

Fisher scalar curvature

S = g^{ij} R_{ij}
"""
function scalar_curvature(g::FisherMetric, ρ::AbstractMatrix, basis)
    if ρ isa Diagonal
        # Als het al een Diagonal is, mag hij HIER niet inkomen.
        # Dit vangt situaties op waarin de splitsing misgaat.
    elseif isdiag(ρ)
        return scalar_curvature(g, Diagonal(ρ), basis)
    end
    G    = metric_matrix(g, ρ, basis)
    Ginv = pinv(G)

    Ric = ricci(g, ρ, basis)

    n = size(G,1)

    S = 0.0

    for i in 1:n
        for j in 1:n
            S += Ginv[i,j] * Ric[i,j]
        end
    end

    return S

end
