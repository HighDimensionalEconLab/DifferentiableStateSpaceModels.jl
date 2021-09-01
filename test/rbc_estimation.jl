using DifferentiableStateSpaceModels,
    LinearAlgebra,
    Optim,
    Test,
    TimerOutputs,
    Turing,
    Zygote,
    BenchmarkTools

Turing.setadbackend(:zygote)
using Turing: @addlogprob!
# Create models from modules and then solve
model_rbc = @include_example_module(Examples.rbc_observables_benchmark)
model_rbc_second = @include_example_module(Examples.rbc_observables_benchmark, 2)

p_f = [0.7, 0.1, 0.01, sqrt(0.0001)]
# Generate artificial data
p = [0.4, 0.96]
sol = generate_perturbation(model_rbc, p; p_f)
sol_second = generate_perturbation(model_rbc_second, p; p_f)
T = 20
ϵ = [randn(model_rbc.n_ϵ) for _ = 1:T]
x0 = zeros(model_rbc.n_x)
fake_z = solve(sol, x0, (0, T), LTI(); noise = ϵ).z
fake_z_second = solve(sol_second, x0, (0, T), QTI(); noise = ϵ).z

## Turing model, Kalman filter
@model function rbc_kalman(z, m, p_f, cache)
    α ~ Uniform(0.2, 0.8)
    β ~ Uniform(0.5, 0.99)
    p = [α, β]
    sol = generate_perturbation(m, p; p_f, cache, settings = PerturbationSolverSettings(;
    print_level = 1))
    if !(sol.retcode == :Success)
        @addlogprob! -Inf
        return
    end
    @addlogprob! solve(sol, sol.x_ergodic, (0, length(z)); observables = z).logpdf
end

turing_model = rbc_kalman(fake_z, model_rbc, p_f, allocate_cache(model_rbc))
n_samples = 20
n_adapts = 5
δ = 0.65
chain = sample(turing_model, NUTS(n_adapts, δ), n_samples; progress = true)
θ_MAP = optimize(turing_model, MAP())

## Turing model, Joint likelihood
@model function rbc_joint(z, m, p_f, cache, x0 = zeros(m.n_x))
    α ~ Uniform(0.2, 0.8)
    β ~ Uniform(0.5, 0.99)
    p = [α, β]
    T = length(z)
    ϵ_draw ~ MvNormal(T, 1.0)
    ϵ = map(i -> ϵ_draw[((i-1)*m.n_ϵ+1):(i*m.n_ϵ)], 1:T)
    # println(p)
    sol = generate_perturbation(m, p; p_f, cache)
    if !(sol.retcode == :Success)
        @addlogprob! -Inf
        return
    end
    @addlogprob! solve(sol, x0, (0, T); noise = ϵ, observables = z).logpdf
end

turing_model = rbc_joint(fake_z, model_rbc, p_f, allocate_cache(model_rbc))
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

## Turing model, Joint likelihood, Second-order
@model function rbc_second(z, m, p_f, cache, x0 = zeros(m.n_x))
    α ~ Uniform(0.2, 0.8)
    β ~ Uniform(0.5, 0.99)
    p = [α, β]
    T = length(z)
    ϵ_draw ~ MvNormal(T, 1.0)
    ϵ = map(i -> ϵ_draw[((i-1)*m.n_ϵ+1):(i*m.n_ϵ)], 1:T)
    # println(p)
    # settings = PerturbationSolverSettings(; print_level = 2)
    sol = generate_perturbation(m, p; p_f, cache)
    if !(sol.retcode == :Success)
        @addlogprob! -Inf
        return
    end
    @addlogprob! solve(sol, x0, (0, T); noise = ϵ, observables = z).logpdf
end
turing_model = rbc_second(fake_z_second, model_rbc_second, p_f, allocate_cache(model_rbc_second))
n_samples = 20
n_adapts = 5
δ = 0.65
max_depth = 2
chain = sample(turing_model, NUTS(n_adapts, δ; max_depth), n_samples; progress = true)
