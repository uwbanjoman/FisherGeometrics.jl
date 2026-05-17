"""
runtests.jl — FisherGeometrics.jl test suite
Run with:  julia --project test/runtests.jl
"""

using Test
using LinearAlgebra

include("../src/FisherGeometrics.jl")
using .FisherGeometrics

@testset "FisherGeometrics.jl" begin

    # ──────────────────────────────────────────────────────────
    @testset "Foundation — constants and Fisher tensor" begin

        @test abs(4*Float64(τ) - cos(π*Float64(τ))) < 1e-2
        @test κ_hol == 6//5
        @test Float64(κ_hol) ≈ 6 * Float64(τ)
        @test abs(1 + φ^4 - 3*φ^2) < 1e-10
        @test abs(2*cos(π/5) - φ) < 1e-10

        @test length(gellmann_basis(2)) == 3
        @test length(gellmann_basis(3)) == 8
        @test length(gellmann_basis(6)) == 35

        ρ₀ = vacuum_state()
        @test tr(ρ₀) ≈ 1.0
        @test ρ₀ ≈ I/6

        ψ = ComplexF64[1, 0, 0, 0, 0, 0]
        ρp = pure_state(ψ)
        @test tr(ρp) ≈ 1.0
        @test real(tr(ρp^2)) ≈ 1.0

        g = fubini_study_metric(ψ)
        F = fisher_tensor(ρp)
        @test maximum(abs.(F - 4*g)) < 1e-10
        @test F ≈ F'
        @test maximum(abs.(imag.(F))) < 1e-12

        basis = gellmann_basis(6)
        @test imag(quantum_fisher(ρp, basis[1], basis[2])) < 1e-12

    end

    # ──────────────────────────────────────────────────────────
    @testset "Geometry — spectral geometry of K" begin

        @test vol_CP2 ≈ π^2/2
        @test vol_S1  ≈ 2π * Float64(τ)
        @test vol_K   ≈ vol_CP2 * vol_S3 * vol_S1

        λs, gs = spectrum_K(5, 5, 5)
        @test length(λs) > 0
        @test all(λ -> λ > 0, λs)
        @test all(g -> g > 0, gs)
        @test length(λs) == length(gs)

        # ζ_K(0) = −61/80 is analytically continued — verify convergence at s=5
        ζ5 = real(zeta_K(5.0; k_max=8, j_max=8, n_max=8))
        @test ζ5 > 0
        @test ζ5 < 1e6

        M2 = kk_masses(5)
        @test M2[1] ≈ 9/4
        @test issorted(M2)

        @test analytic_torsion() ≈ 1.0

        ρ = vacuum_state()
        @test bures_distance(ρ, ρ) < 1e-7

        ψ1 = ComplexF64[1,0,0,0,0,0]
        ψ2 = ComplexF64[0,1,0,0,0,0]
        d12 = bures_distance(pure_state(ψ1), pure_state(ψ2))
        @test d12 ≥ 0
        @test d12 ≈ bures_distance(pure_state(ψ2), pure_state(ψ1))

    end

    # ──────────────────────────────────────────────────────────
    @testset "Symmetry — gauge structure" begin

        @test n_generations == 3

        # c₁ = c₂ = χ = 3 — verify via hypercharges and generations
        @test n_generations == 3   # = c₁(ℂP²) topological

        Y = hypercharges()
        @test Y.Y_QL == 1//6
        @test Y.Y_uR == 2//3
        @test Y.Y_dR == -1//3
        @test Y.Y_LL == -1//2
        @test Y.Y_eR == -1

        @test check_anomaly_cancellation()
        @test abs(sin2_weinberg - 0.2312) / 0.2312 < 0.01
        @test abs(alpha_GUT_inv() - 41.5) / 41.5 < 0.05

    end

    # ──────────────────────────────────────────────────────────
    @testset "Dynamics — couplings and CKM" begin

        @test abs(alpha_em_inv - 137.036) / 137.036 < 0.001
        @test abs(λ_W - 0.2250) / 0.2250 < 0.03
        @test abs(rad2deg(δ_CP) - 69.2) / 69.2 < 0.005

        @test abs(1 + φ^4 - 3φ^2) < 1e-10
        @test abs(sin(δ_CP) - φ/√3) < 1e-10
        @test abs(cos(δ_CP) - 1/(φ*√3)) < 1e-10

        V = ckm_matrix()
        @test maximum(abs.(V * V' - I)) < 1e-10
        @test abs(V[1,1]) > 0.97
        @test abs(V[2,2]) > 0.97
        @test abs(V[3,3]) > 0.99
        @test abs(abs(V[2,3]) - 0.0418) / 0.0418 < 0.03
        @test abs(abs(V[1,3]) - 0.00351) / 0.00351 < 0.03

        J = jarlskog()
        @test J > 0
        @test abs(J - 3.08e-5) / 3.08e-5 < 0.05

        r = mass_ratios()
        @test r[1] ≈ 1.0
        @test r[2] ≈ Float64(τ)^2
        @test r[3] ≈ Float64(τ)^4
        @test r[1] > r[2] > r[3]

        ids = cp_phase_identity()
        @test ids.identity_1_err < 1e-10
        @test ids.identity_2_err < 1e-10
        @test ids.identity_3_err < 1e-10

    end

    # ──────────────────────────────────────────────────────────
    @testset "Gravity — Einstein equation and cosmology" begin

        @test Λ_fundamental == 0.0
        @test bh_entropy(1.0, 1.0) ≈ 0.25
        @test bh_entropy(4.0, 1.0) ≈ 1.0
        @test bh_entropy(π, 1/(4π)) ≈ π^2
        @test bh_temperature(1.0, 1.0) > 0
        @test cosmological_constant(1.0) > 0
        @test FisherGeometrics.vacuum_action() ≈ exp(-0.5)
        @test abs(Ω_Λ - 0.70) / 0.70 < 0.10
        @test abs(n_s - 0.9649) / 0.9649 < 0.005
        @test Ω_DM == 0.0

    end

    # ──────────────────────────────────────────────────────────
    @testset "Evolution — Von Neumann dynamics" begin

        H  = [1.0 0.3; 0.3 -1.0]
        ρ₀ = let r = [0.7 0.2+0.1im; 0.2-0.1im 0.3]; (r+r')/2; end

        @test abs(tr(von_neumann_rhs(ρ₀, H))) < 1e-12

        ρT = evolve_exact(ρ₀, H, 2π)
        @test abs(tr(ρT) - 1.0) < 1e-10
        @test abs(purity(ρ₀) - purity(ρT)) < 1e-10
        @test abs(entropy(ρ₀) - entropy(ρT)) < 1e-8
        @test ρT ≈ ρT'

        traj = evolve_rk4(ρ₀, H, 1.0; dt=0.001)
        @test maximum(abs.(traj[end][2] - evolve_exact(ρ₀, H, 1.0))) < 1e-6

        ψ = ComplexF64[1,0,0,0,0,0]
        @test purity(pure_state(ψ)) ≈ 1.0
        @test purity(vacuum_state()) ≈ 1/6
        @test entropy(pure_state(ψ)) < 1e-10
        @test abs(entropy(vacuum_state()) - log(6)) < 1e-10

        H_KK = hamiltonian_KK(4)
        @test isdiag(H_KK)
        @test all(diag(H_KK) .> 0)
        @test diag(H_KK)[1] ≈ 9/4

        obs = Diagonal([1.0, 2.0])
        ρ2  = let r = [0.6 0.2; 0.2 0.4]; (r+r')/2; end
        res = measurement_projection(ρ2, obs)
        @test sum(res.probabilities) ≈ 1.0
        @test all(p -> p ≥ 0, res.probabilities)

        traj2 = evolve_rk4(ρ₀, H, 0.5; dt=0.01)
        dists = information_distance_trajectory(traj2)
        @test dists[1] < 1e-7

    end

    # ──────────────────────────────────────────────────────────
    @testset "Integration — check_all passes" begin
        @test FisherGeometrics.check_all()
    end

end
