
using LinearAlgebra
using FisherGeometrics

n_val = 6
T_basis = su_basis(n_val)
N = length(T_basis)

# Metric at ρ* (diagonal: g_{ab} = (n/8) δ_{ab})
g_num   = (n_val/8) * Matrix{Float64}(I, N, N)
g_inv   = (8/n_val) * Matrix{Float64}(I, N, N)

# Bures vacuum metric g_{ab} = (n/8)δ_{ab} at ρ* = I/6
g   = (6/8) * Matrix(I, 35, 35)

# Stap by stap
d   = d_tensor(T_basis)
Γ   = christoffel_vacuum(d, 6)
Q   = riemann_quadratic(Γ)
Ric = ricci_tensor(Q)
S_F = ricci_scalar(Ric, g_inv)
G   = einstein_tensor(Ric, g, S_F)

# Or all in one go
FG = bures_einstein(T_basis)
FG.S_F
FG.λ
FG.Λ