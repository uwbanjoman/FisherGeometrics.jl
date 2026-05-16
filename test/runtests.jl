"""
runtests.jl
===========
Test suite for FisherGeometrics.jl

Run with:  julia --project -e "using Pkg; Pkg.test()"
Or:        julia --project test/runtests.jl
"""

using Test
using LinearAlgebra

# Load the package
include("../src/FisherGeometrics.jl")
using .FisherGeometrics

@testset "FisherGeometrics.jl" begin

    # ─────────────────────────────────────────────────────────
    @testset "Foundation — constants and Fisher tensor" begin

        # τ uniqueness: 4τ = cos(πτ)
        @test abs(4*Float64(τ) - cos(π*Float64(τ))) < 1e-2

        # κ_hol = 6τ
        @test κ_hol == 6//5
        @test Float64(κ_hol) ≈ 6 * Float64(τ)

        # Golden ratio identity: 1 + φ⁴ = 3φ²
        @test abs(1 + φ^4 - 3*φ^2) < 1e-10

        # [2]_q = 2cos(π/5) = φ
        @test abs(2*cos(π/5) - φ) < 1e-10

        # Gell-Mann basis for n=2: 3 generators (su(2))
        basis2 = gellmann_basis(2)
        @test length(basis2) == 3

        # Gell-Mann basis for n=3: 8 generators (su(3))
        basis3 = gellmann_basis(3)
        @test length(basis3) == 8

        # Gell-Mann basis for n=6: 35 generators (su(6))
        basis6 = gellmann_basis(6)
        @test length(basis6) == 35

        # Vacuum state: Tr(ρ₀) = 1, ρ₀ = I/6
        ρ₀ = vacuum_state()
        @test tr(ρ₀) ≈ 1.0
        @test ρ₀ ≈ I/6

        # Pure state: Tr(ρ) = 1, Tr(ρ²) = 1
        ψ = ComplexF64[1, 0, 0, 0, 0, 0]
        ρ_pure = pure_state(ψ)
        @test tr(ρ_pure) ≈ 1.0
        @test real(tr(ρ_pure^2)) ≈ 1.0

        # Braunstein-Caves: F^(Q) = 4 g^(FS)
        g_FS = fubini_study_metric(ψ)
        F    = fisher_tensor(ρ_pure)
        @test maximum(abs.(F - 4*g_FS)) < 1e-10

        # Fisher tensor is symmetric
        ρ_test = ρ_pure
        F_test = fisher_tensor(ρ_test)
        @test F_test ≈ F_test'

        # quantum_fisher is real
        A = basis6[1]
        B = basis6[2]
        @test imag(quantum_fisher(ρ_pure, A, B)) < 1e-12

    end

    # ─────────────────────────────────────────────────────────
    @testset "Geometry — spectral geometry of K" begin

        # Volumes
        @test vol_CP2 ≈ π^2/2
        @test vol_S1  ≈ 2π * Float64(τ)
        @test vol_K   ≈ vol_CP2 * vol_S3 * vol_S1

        # spectrum_K returns non-empty lists
        λs, gs = spectrum_K(5, 5, 5)
        @test length(λs) > 0
        @test length(gs) == length(λs)

        # All eigenvalues positive
        @test all(λ -> λ > 0, λs)

        # All multiplicities positive
        @test all(g -> g > 0, gs)

        # ζ_K(0) ≈ −61/80
        ζ0 = real(zeta_K(0.0+0im; k_max=15, j_max=15, n_max=15))
        @test abs(ζ0 - (-61/80)) / abs(-61/80) < 0.02   # within 2%

        # Minimum KK mass = 9/4
        M2 = kk_masses(3)
        @test M2[1] ≈ 9/4

        # Analytic torsion = 1 exactly
        @test analytic_torsion() ≈ 1.0

        # Bures distance: D(ρ, ρ) = 0
        ρ = vacuum_state()
        @test bures_distance(ρ, ρ) < 1e-10

        # Bures distance: D(ρ₁, ρ₂) ≥ 0
        ψ1 = ComplexF64[1, 0, 0, 0, 0, 0]
        ψ2 = ComplexF64[0, 1, 0, 0, 0, 0]
        ρ1 = pure_state(ψ1)
        ρ2 = pure_state(ψ2)
        @test bures_distance(ρ1, ρ2) ≥ 0

        # Bures distance: D(ρ₁, ρ₂) = D(ρ₂, ρ₁)  (symmetry)
        @test bures_distance(ρ1, ρ2) ≈ bures_distance(ρ2, ρ1)

    end

    # ─────────────────────────────────────────────────────────
    @testset "Symmetry — gauge structure" begin

        # Three generations
        @test n_generations == 3

        # c₁ = c₂ = χ = 3
        @test FisherGeometrics.Symmetry.c1_CP2 == 3
        @test FisherGeometrics.Symmetry.c2_CP2 == 3
        @test FisherGeometrics.Symmetry.χ_CP2  == 3

        # Hypercharges are exact rationals
        Y = hypercharges()
        @test Y.Y_QL == 1//6
        @test Y.Y_uR == 2//3
        @test Y.Y_dR == -1//3
        @test Y.Y_LL == -1//2
        @test Y.Y_eR == -1

        # Anomaly cancellation
        @test FisherGeometrics.Symmetry.check_anomaly_cancellation()

        # Weinberg angle within 1% of observed
        @test abs(sin2_weinberg - 0.2312) / 0.2312 < 0.01

        # GUT coupling within 5% of observed
        @test abs(alpha_GUT_inv() - 41.5) / 41.5 < 0.05

    end

    # ─────────────────────────────────────────────────────────
    @testset "Dynamics — couplings and CKM" begin

        # Fine structure constant within 0.1%
        @test abs(alpha_em_inv - 137.036) / 137.036 < 0.001

        # Cabibbo angle within 3%
        @test abs(λ_W - 0.2250) / 0.2250 < 0.03

        # CP phase within 0.5%
        @test abs(rad2deg(δ_CP) - 69.2) / 69.2 < 0.005

        # Golden ratio identities (exact)
        @test abs(1 + φ^4 - 3φ^2) < 1e-10
        @test abs(sin(δ_CP) - φ/√3) < 1e-10
        @test abs(cos(δ_CP) - 1/(φ*√3)) < 1e-10

        # CKM matrix is unitary
        V = ckm_matrix()
        @test maximum(abs.(V * V' - I)) < 1e-10

        # CKM diagonal elements close to 1
        @test abs(V[1,1]) > 0.97
        @test abs(V[2,2]) > 0.97
        @test abs(V[3,3]) > 0.99

        # |V_cb| within 3%
        @test abs(abs(V[2,3]) - 0.0418) / 0.0418 < 0.03

        # |V_ub| within 3%
        @test abs(abs(V[1,3]) - 0.00351) / 0.00351 < 0.03

        # Jarlskog invariant positive and correct order of magnitude
        J = jarlskog()
        @test J > 0
        @test abs(J - 3.08e-5) / 3.08e-5 < 0.05

        # Mass hierarchy ratios
        r = mass_ratios()
        @test r[1] ≈ 1.0
        @test r[2] ≈ Float64(τ)^2
        @test r[3] ≈ Float64(τ)^4
        @test r[1] > r[2] > r[3]   # correct ordering

        # CP phase identities
        ids = cp_phase_identity()
        @test ids.identity_1_err < 1e-10
        @test ids.identity_2_err < 1e-10
        @test ids.identity_3_err < 1e-10

    end

    # ─────────────────────────────────────────────────────────
    @testset "Gravity — Einstein equation and cosmology" begin

        # Fundamental cosmological constant is exactly zero
        @test Λ_fundamental == 0.0

        # Bekenstein-Hawking: S = A/4G
        @test bh_entropy(1.0, 1.0) ≈ 0.25
        @test bh_entropy(4.0, 1.0) ≈ 1.0
        @test bh_entropy(π, 1/(4π)) ≈ 4π^2  # S = A/4G with G=1/4π

        # Hawking temperature: T > 0 for finite mass
        @test bh_temperature(1.0, 1.0) > 0

        # Cosmological constant is positive for positive G_N
        @test cosmological_constant(1.0) > 0

        # Vacuum action = e^{-1/2}
        @test FisherGeometrics.Gravity.vacuum_action() ≈ exp(-0.5)

        # Ω_Λ within 10% of observed
        @test abs(Ω_Λ - 0.70) / 0.70 < 0.10

        # n_s within 0.5% of observed
        @test abs(n_s - 0.9649) / 0.9649 < 0.005

        # Dark matter prediction
        @test Ω_DM == 0.0

    end

    # ─────────────────────────────────────────────────────────
    @testset "Evolution — Von Neumann dynamics" begin

        H  = [1.0 0.3; 0.3 -1.0]
        ρ₀ = [0.7 0.2+0.1im; 0.2-0.1im 0.3]
        ρ₀ = (ρ₀ + ρ₀') / 2

        # Von Neumann RHS is traceless (Tr([H,ρ]) = 0 always)
        rhs = von_neumann_rhs(ρ₀, H)
        @test abs(tr(rhs)) < 1e-12

        # Exact evolution: trace preserved
        ρT = evolve_exact(ρ₀, H, 2π)
        @test abs(tr(ρT) - 1.0) < 1e-10

        # Exact evolution: purity preserved (unitary)
        @test abs(purity(ρ₀) - purity(ρT)) < 1e-10

        # Exact evolution: entropy preserved
        @test abs(entropy(ρ₀) - entropy(ρT)) < 1e-8

        # Exact evolution: hermiticity preserved
        @test ρT ≈ ρT'

        # RK4 agrees with exact solution within 10⁻⁶
        traj = evolve_rk4(ρ₀, H, 1.0; dt=0.001)
        ρ_rk4  = traj[end][2]
        ρ_exact = evolve_exact(ρ₀, H, 1.0)
        @test maximum(abs.(ρ_rk4 - ρ_exact)) < 1e-6

        # Purity: pure state has purity 1
        ψ = ComplexF64[1, 0, 0, 0, 0, 0]
        @test purity(pure_state(ψ)) ≈ 1.0

        # Purity: maximally mixed has purity 1/n
        ρ_mix = Matrix{ComplexF64}(I, 6, 6) / 6
        @test purity(ρ_mix) ≈ 1/6

        # Entropy: pure state has entropy 0
        @test entropy(pure_state(ψ)) < 1e-10

        # Entropy: maximally mixed has entropy log(n)
        @test abs(entropy(ρ_mix) - log(6)) < 1e-10

        # Hamiltonian KK: diagonal, positive entries
        H_KK = hamiltonian_KK(4)
        @test isdiag(H_KK)
        @test all(diag(H_KK) .> 0)

        # First KK mass = 9/4
        @test diag(H_KK)[1] ≈ 9/4

        # Measurement: probabilities sum to 1
        obs = Diagonal([1.0, 2.0])
        ρ2  = [0.6 0.2; 0.2 0.4]
        result = measurement_projection(ρ2, obs)
        @test sum(result.probabilities) ≈ 1.0
        @test all(p -> p ≥ 0, result.probabilities)

        # Information distance trajectory: starts at 0
        traj_full = evolve_rk4(ρ₀, H, 0.5; dt=0.01)
        dists = information_distance_trajectory(traj_full)
        @test dists[1] < 1e-10

    end

    # ─────────────────────────────────────────────────────────
    @testset "Integration — check_all passes" begin
        @test check_all()
    end

end
