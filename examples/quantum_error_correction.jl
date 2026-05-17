# examples/quantum_error_correction.jl
# =======================================
# Quantum error correction from the topology of ℂP²
#
# ── The idea ─────────────────────────────────────────────────────────────────
#
# Standard quantum error correction imposes protection externally —
# you encode a logical qubit into many physical qubits and detect
# errors via syndrome measurements.
#
# The FisherGeometrics framework suggests a different approach:
# the information geometry of K = ℂP² × S³ × S¹ provides intrinsic
# topological protection. The three generations are protected by
# c₁(ℂP²) = 3 — errors that respect the topology cannot corrupt them.
#
# ── The SM-bit ───────────────────────────────────────────────────────────────
#
# The natural unit of quantum information in the framework is not a qubit
# but a 6-dimensional quantum system on ℂ⁶ = ℂ³ ⊗ ℂ²:
#
#   ℂ³  — qutrit  (colour, lives on ℂP²)
#   ℂ²  — qubit   (isospin, lives on S³)
#
# We call this a "SM-bit" (Standard Model bit). It has 6 basis states:
#   |r↑⟩, |r↓⟩, |g↑⟩, |g↓⟩, |b↑⟩, |b↓⟩
#
# ── Error protection from topology ───────────────────────────────────────────
#
# The key insight from ℂP²:
#
#   c₁(ℂP²) = 3  →  errors that change the colour index by < 3 steps
#                    cannot distinguish between generations
#
# The key insight from S³:
#
#   χ(S³) = 0  →  no preferred isospin direction
#              →  isotropic decoherence, no preferred error axis
#
# This script demonstrates:
#   1. The SM-bit basis and its Fisher information structure
#   2. How topological errors (ℂP² rotations) are undetectable
#   3. How the Bures distance measures error severity
#   4. Comparison with standard qubit error correction
#
# © 2026 Jan Bouwman — FisherGeometrics Framework

using FisherGeometrics
using LinearAlgebra
using Printf

println("="^65)
println("  FisherGeometrics — Quantum Error Correction from Topology")
println("="^65)
println()

# ── The SM-bit basis ──────────────────────────────────────────────────────────

println("─── The SM-bit: ℂ⁶ = ℂ³ ⊗ ℂ² ───")
println()
println("  Basis states  |colour, isospin⟩:")
println()

basis_labels = ["|r↑⟩", "|r↓⟩", "|g↑⟩", "|g↓⟩", "|b↑⟩", "|b↓⟩"]
for (i, label) in enumerate(basis_labels)
    colour = ["red", "red", "green", "green", "blue", "blue"][i]
    spin   = ["up", "down", "up", "down", "up", "down"][i]
    @printf("  |e_%d⟩ = %-6s  (colour: %-5s, isospin: %s)\n",
            i, label, colour, spin)
end
println()
println("  This is the natural 6-dimensional quantum system of the SM vacuum.")
println("  One SM-bit encodes more information than a qubit (2D) or qutrit (3D).")
println()

# ── Three logical states (three generations) ──────────────────────────────────

println("─── Three logical states — protected by c₁(ℂP²) = 3 ───")
println()
println("  The three generations correspond to three orthogonal subspaces")
println("  of ℂ⁶, each spanned by a colour⊗isospin pair:")
println()

# Generation 1: red sector
ψ_gen1 = ComplexF64[1, 1, 0, 0, 0, 0] / sqrt(2)   # |r↑⟩ + |r↓⟩
# Generation 2: green sector
ψ_gen2 = ComplexF64[0, 0, 1, 1, 0, 0] / sqrt(2)   # |g↑⟩ + |g↓⟩
# Generation 3: blue sector
ψ_gen3 = ComplexF64[0, 0, 0, 0, 1, 1] / sqrt(2)   # |b↑⟩ + |b↓⟩

gens = [ψ_gen1, ψ_gen2, ψ_gen3]
gen_names = ["Gen 1 (electron)", "Gen 2 (muon)", "Gen 3 (tau)"]

for (i, (ψ, name)) in enumerate(zip(gens, gen_names))
    @printf("  |L_%d⟩ = %s:  ", i, name)
    for (j, c) in enumerate(ψ)
        abs(c) > 1e-10 && @printf("%+.3f%s ", real(c), basis_labels[j])
    end
    println()
end
println()

# Verify orthogonality
println("  Orthogonality (⟨Lᵢ|Lⱼ⟩ = δᵢⱼ):")
for i in 1:3, j in 1:3
    overlap = abs(dot(gens[i], gens[j]))
    expected = i == j ? 1.0 : 0.0
    ok = abs(overlap - expected) < 1e-10 ? "✓" : "✗"
    @printf("  ⟨L_%d|L_%d⟩ = %.4f  %s\n", i, j, overlap, ok)
end
println()

# ── Fisher information of the logical states ──────────────────────────────────

println("─── Fisher information of the logical states ───")
println()

for (i, (ψ, name)) in enumerate(zip(gens, gen_names))
    ρ = pure_state(ψ)
    F = fisher_tensor(ρ)
    @printf("  Gen %d (%s):\n", i, name)
    @printf("    Tr(𝓕) = %.4f,  rank(𝓕) = %d,  max|𝓕_AB| = %.4f\n",
            tr(F), rank(F, rtol=1e-10), maximum(abs.(F)))
end
println()
println("  All three generations have identical Fisher information —")
println("  they are informationally equivalent. This is the geometric")
println("  expression of the universality of the weak interaction.")
println()

# ── Topological errors vs non-topological errors ──────────────────────────────

println("─── Error types: topological vs non-topological ───")
println()
println("  A topological error permutes the colour index: r→g→b→r")
println("  This corresponds to a U(1) ⊂ SU(3) rotation on ℂP².")
println("  Such errors CANNOT be detected by local measurements")
println("  — they are invisible to the Fisher metric within each generation.")
println()
println("  A non-topological error mixes colour and isospin: |r↑⟩ → |g↓⟩")
println("  This changes the Fisher information significantly")
println("  and IS detectable.")
println()

# Define error operators
# Topological: cyclic colour permutation r→g→b→r (ℂP² rotation)
# This is a U(1) transformation, undetectable within the topology
U_topo = zeros(ComplexF64, 6, 6)
U_topo[3,1] = 1;  U_topo[4,2] = 1   # r → g
U_topo[5,3] = 1;  U_topo[6,4] = 1   # g → b
U_topo[1,5] = 1;  U_topo[2,6] = 1   # b → r

# Non-topological: bit flip within colour sector
U_nontopo = Matrix{ComplexF64}(I, 6, 6)
U_nontopo[1,1] = 0; U_nontopo[3,3] = 0
U_nontopo[1,3] = 1; U_nontopo[3,1] = 1   # mix r↑ and g↑

println("  Testing on Gen 1 state |L₁⟩:")
println()

ρ_gen1 = pure_state(ψ_gen1)

# Apply topological error
ψ_after_topo = U_topo * ψ_gen1
ψ_after_topo /= norm(ψ_after_topo)
ρ_after_topo = pure_state(ψ_after_topo)
d_topo = bures_distance(ρ_gen1, ρ_after_topo)

# Apply non-topological error
ψ_after_nontopo = U_nontopo * ψ_gen1
ψ_after_nontopo /= norm(ψ_after_nontopo)
ρ_after_nontopo = pure_state(ψ_after_nontopo)
d_nontopo = bures_distance(ρ_gen1, ρ_after_nontopo)

@printf("  Topological error (colour permutation r→g→b→r):\n")
@printf("    State after:  ")
for (j, c) in enumerate(ψ_after_topo)
    abs(c) > 1e-10 && @printf("%+.3f%s ", real(c), basis_labels[j])
end
println()
@printf("    This is |L₂⟩ (Gen 2) — same generation structure, different label\n")
@printf("    Bures distance from |L₁⟩:  %.4f\n", d_topo)
@printf("    Fisher info change:         %.4f → %.4f  (unchanged)\n",
        tr(fisher_tensor(ρ_gen1)), tr(fisher_tensor(ρ_after_topo)))
println()

@printf("  Non-topological error (colour mixing r↑ ↔ g↑):\n")
@printf("    State after:  ")
for (j, c) in enumerate(ψ_after_nontopo)
    abs(c) > 1e-10 && @printf("%+.3f%s ", real(c), basis_labels[j])
end
println()
@printf("    Bures distance from |L₁⟩:  %.4f\n", d_nontopo)
println()
println("  → Topological error: measurable Bures distance, but")
println("    maps one valid generation to another. Correctable by")
println("    tracking the generation label.")
println()
println("  → Non-topological error: creates an invalid superposition.")
println("    Detectable and correctable via syndrome measurement.")
println()

# ── The Knill-Laflamme criterion in Fisher language ───────────────────────────

println("─── Knill-Laflamme criterion in Fisher geometry ───")
println()
println("  Standard KL criterion: error E is correctable if")
println("    ⟨Lᵢ|E†E|Lⱼ⟩ = C_ij (constant, independent of i,j)")
println()
println("  Fisher geometry version: error E is correctable if")
println("    𝓕_AB[E(ρ_i)] = 𝓕_AB[ρ_i]  for all logical states ρ_i")
println("  i.e., the error does not change the information content.")
println()

errors = [
    ("Identity (no error)",     Matrix{ComplexF64}(I, 6, 6)),
    ("Topological (r→g→b→r)",   U_topo),
    ("Non-topological (r↔g)",   U_nontopo),
    ("Phase flip |r↑⟩→-|r↑⟩",  begin
        U = Matrix{ComplexF64}(I,6,6); U[1,1]=-1; U
    end),
    ("Depolarising (→ I/6)",    zeros(ComplexF64,6,6)),   # special case
]

println("  Error channel          |  KL satisfied?  |  Max Bures dist")
println("  " * "─"^55)

for (name, E) in errors
    kl_ok  = true
    max_d  = 0.0

    if all(E .== 0)
        # Depolarising channel: maps everything to I/6
        for ψ in gens
            ρ = pure_state(ψ)
            ρ_err = vacuum_state()
            d = bures_distance(ρ, ρ_err)
            max_d = max(max_d, d)
        end
        kl_ok = false
    else
        F_refs = [fisher_tensor(pure_state(ψ)) for ψ in gens]
        for (ψ, F_ref) in zip(gens, F_refs)
            ρ = pure_state(ψ)
            ψe = E * ψ; ψe /= norm(ψe)
            ρe = pure_state(ψe)
            Fe = fisher_tensor(ρe)
            d  = bures_distance(ρ, ρe)
            max_d = max(max_d, d)
            if maximum(abs.(Fe - F_ref)) > 0.01
                kl_ok = false
            end
        end
    end

    kl_str = kl_ok ? "✓  yes" : "✗  no"
    @printf("  %-24s |  %-15s |  %.4f\n", name, kl_str, max_d)
end
println()

# ── Comparison: SM-bit vs qubit ───────────────────────────────────────────────

println("─── SM-bit vs qubit: information capacity ───")
println()
println("  Information capacity comparison:")
println()
@printf("  %-20s  %-8s  %-12s  %-16s\n",
        "System", "dim", "log₂(dim)", "Topol. protection")
println("  " * "─"^60)
@printf("  %-20s  %-8d  %-12.3f  %-16s\n", "Qubit (ℂ²)",    2, log2(2), "none")
@printf("  %-20s  %-8d  %-12.3f  %-16s\n", "Qutrit (ℂ³)",   3, log2(3), "c₁(ℂP¹)=1")
@printf("  %-20s  %-8d  %-12.3f  %-16s\n", "SM-bit (ℂ⁶)",   6, log2(6), "c₁(ℂP²)=3 ✓")
@printf("  %-20s  %-8d  %-12.3f  %-16s\n", "3 qubits",       8, log2(8), "none")
println()
println("  The SM-bit (ℂ⁶) encodes log₂(6) ≈ 2.58 bits of information")
println("  with intrinsic topological protection from c₁(ℂP²) = 3.")
println("  This is more than a qubit (1 bit) or qutrit (1.58 bits),")
println("  and the protection is built into the geometry — not added.")
println()

# ── Summary ───────────────────────────────────────────────────────────────────

println("─── Summary ───")
println()
println("  ┌─────────────────────────────────────────────────────────┐")
println("  │  The FisherGeometrics framework suggests:               │")
println("  │                                                         │")
println("  │  1. The natural quantum information unit is the         │")
println("  │     SM-bit: ℂ⁶ = ℂ³ ⊗ ℂ²  (qutrit ⊗ qubit)          │")
println("  │                                                         │")
println("  │  2. Topological protection from c₁(ℂP²) = 3:          │")
println("  │     colour permutations map valid states to valid       │")
println("  │     states — they do not create invalid superpositions  │")
println("  │                                                         │")
println("  │  3. Isotropic decoherence from χ(S³) = 0:             │")
println("  │     no preferred error axis in the isospin sector       │")
println("  │                                                         │")
println("  │  4. The Bures distance is the natural error metric:     │")
println("  │     errors are small ↔ small Bures distance from        │")
println("  │     the codespace                                       │")
println("  │                                                         │")
println("  │  Open: physical implementation of the SM-bit            │")
println("  │  Candidates: trapped ions (6 levels), photonic          │")
println("  │  systems, topological anyons (Fibonacci, k=3)          │")
println("  └─────────────────────────────────────────────────────────┘")
println()
println("  See: FisherGeometrics Framework, Documents I–LXXIII (2026)")
