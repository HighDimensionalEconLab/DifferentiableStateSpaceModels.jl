using DifferentiableStateSpaceModels, LinearAlgebra, Optim, Test, Turing, Zygote
using DifferentiableStateSpaceModels.Examples
using Turing: @addlogprob!

Turing.setadbackend(:zygote)

# Create models from modules and then solve
# p_d = [α, β]
model_rbc = @include_example_module(Examples.rbc_observables)

p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.01)
p_d = (α=0.5, β=0.95)
p_d_input = [0.5, 0.95]
sol = generate_perturbation(model_rbc, p_d, p_f, Val(1))
sol_second = generate_perturbation(model_rbc, p_d, p_f, Val(2))

T = 20
ϵ = [randn(model_rbc.n_ϵ) for _ = 1:T]
x0 = zeros(model_rbc.n_x)
fake_z = solve(sol, x0, (0, T), DifferentiableStateSpaceModels.LTI(); noise = ϵ).z
fake_z_second = solve(sol_second, x0, (0, T), DifferentiableStateSpaceModels.QTI(); noise = ϵ).z

# Turing model, Kalman filter
@model function rbc_kalman(z, m, p_f, cache)
    α ~ Uniform(0.2, 0.8)
    β ~ Uniform(0.5, 0.99)
    p_d = (α = α, β = β)
    # println(p_d)
    sol = generate_perturbation(m, p_d, p_f, Val(1); cache)
    if !(sol.retcode == :Success)
        @addlogprob! -Inf
        return
    end
    @addlogprob! solve(sol, sol.x_ergodic, (0, length(z)); observables = z).logpdf
end

c = SolverCache(model_rbc, Val(1), p_d)
turing_model = rbc_kalman(fake_z, model_rbc, p_f, c)
n_samples = 20
n_adapts = 5
δ = 0.65
chain = sample(turing_model, NUTS(n_adapts, δ), n_samples; progress = true)
θ_MAP = optimize(turing_model, MAP())

# Turing model, Joint likelihood
@model function rbc_joint(z, m, p_f, cache, x0 = zeros(m.n_x))
    α ~ Uniform(0.2, 0.8)
    β ~ Uniform(0.5, 0.99)
    p_d = (α = α, β = β)
    T = length(z)
    ϵ_draw ~ MvNormal(T, 1.0)
    ϵ = map(i -> ϵ_draw[((i-1)*m.n_ϵ+1):(i*m.n_ϵ)], 1:T)
    # println(p_d)
    sol = generate_perturbation(m, p_d, p_f, Val(1); cache)
    if !(sol.retcode == :Success)
        @addlogprob! -Inf
        return
    end
    @addlogprob! solve(sol, x0, (0, T); noise = ϵ, observables = z).logpdf
end

c = SolverCache(model_rbc, Val(1), p_d)
turing_model = rbc_joint(fake_z, model_rbc, p_f, c)
n_samples = 20
n_adapts = 3
δ = 0.65
chain = sample(turing_model, NUTS(n_adapts, δ), n_samples; progress = true)

ϵ_leapfrog = 0.02
n_depth = 2
chain = sample(
    turing_model,
    Gibbs(HMC(ϵ_leapfrog, n_depth, :α, :β), HMC(ϵ_leapfrog, n_depth, :ϵ_draw)),
    n_samples;
    progress = true
)

# Gradient check codes
#=
density = Turing.OptimLogDensity(turing_model, DynamicPPL.DefaultContext())
g = gradient(density, θ_MAP.values.array)
g = gradient(density, [p;vcat(ϵ...)])
=#
# FiniteDiff.finite_difference_gradient(density, [p;vcat(ϵ...)])

# Turing model, Joint likelihood, Second-order
@model function rbc_second(z, m, p_f, cache, x0 = zeros(m.n_x))
    α ~ Uniform(0.2, 0.8)
    β ~ Uniform(0.5, 0.99)
    p_d = (α = α, β = β)
    T = length(z)
    ϵ_draw ~ MvNormal(T, 1.0)
    ϵ = map(i -> ϵ_draw[((i-1)*m.n_ϵ+1):(i*m.n_ϵ)], 1:T)
    # println(p_d)
    sol = generate_perturbation(m, p_d, p_f, Val(2); cache)
    if !(sol.retcode == :Success)
        @addlogprob! -Inf
        return
    end
    @addlogprob! solve(sol, x0, (0, T); noise = ϵ, observables = z).logpdf
end
c = SolverCache(model_rbc, Val(2), p_d)
turing_model = rbc_joint(fake_z, model_rbc, p_f, c)
n_samples = 20
n_adapts = 5
δ = 0.65
max_depth = 2
chain = sample(turing_model, NUTS(n_adapts, δ; max_depth), n_samples; progress = true)
