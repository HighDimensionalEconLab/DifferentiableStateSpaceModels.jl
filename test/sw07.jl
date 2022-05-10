using DifferentiableStateSpaceModels, SparseArrays, LinearAlgebra, Parameters, Test,
      TimerOutputs, BenchmarkTools

isdefined(Main, :SW07) || include(joinpath(pkgdir(DifferentiableStateSpaceModels), "test/generated_models/SW07.jl"))
const m_sw = PerturbationModel(Main.SW07)

# @testset "Dense SW First Order" begin
#     m = @include_example_module(Examples.SW07)

#     p = [10, 0.51, 10, 0, 0.7, 0.742, 0, 0, 0.24, 0.2696, 6.0144, 0.025, 1.5, 0.6361, 1.5,
#          0.3243, 0.8087, 0.47, 0.6, 1.9423, 1.5, 1.488, 0.2347, 0.0593, 0.8762, 0.9977,
#          0.5799, 0.9957, 0.7165, 0, 0, 0, 0.3982, 0.18]
#     generate_perturbation(m, p)
#     reset_timer!()
#     sol = generate_perturbation(m, p)
#     print_timer()
#     @test sol.retcode == :Success
# end

p_d = (ε_w = 10, ρ_ga = 0.51, ε_p = 10, l_bar = 0, Π_bar = 0.7, B = 0.742, μ_w = 0, μ_p = 0, α = 0.24, ψ = 0.2696, φ = 6.0144, δ = 0.025, σ_c = 1.5, λ = 0.6361, ϕ_p = 1.5, ι_w = 0.3243, ξ_w = 0.8087, ι_p = 0.47, ξ_p = 0.6, σ_l = 1.9423, ϕ_w = 1.5, r_π = 1.488, r_Δy = 0.2347, r_y = 0.0593, ρ = 0.8762, ρ_a = 0.9977, ρ_b = 0.5799, ρ_g = 0.9957, ρ_i = 0.7165, ρ_r = 0, ρ_p = 0, ρ_w = 0, γ_bar = 0.3982, gy_ss = 0.18, se_a = 0.4618, se_b = 1.8513, se_g = 0.6090, se_i = 0.6017, se_m = 0.2397, se_π = 0.1455, se_w = 0.2089)
p_f = (Ω_ii = sqrt(1e-5),)

c = SolverCache(m_sw, Val(1), p_d)
sol = generate_perturbation(m_sw, p_d, p_f; cache = c)
