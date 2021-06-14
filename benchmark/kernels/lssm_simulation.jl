using BenchmarkTools, Zygote, LinearAlgebra, Parameters, RecursiveArrayTools
using Zygote: @adjoint

struct ExampleLTI{T1,T2,T3}
    A::T1
    B::T2
    C::T3
end
@adjoint function ExampleLTI(A, B, C)
    return ExampleLTI(A, B, C), Δ -> (Δ.A, Δ.B, Δ.C)
end

abstract type AbstractStateSpaceAlgorithm end
struct GeneralStateSpace <: AbstractStateSpaceAlgorithm end
struct LTIStateSpace <: AbstractStateSpaceAlgorithm end

# General algorithm (nonlinear, time-varying)
function solve(f, g, h, p, u0, noise, tspan, alg::GeneralStateSpace)
    T = tspan[2]
    u1 = f(p, u0, 1) .+ g(p, u0, 1) * noise[1]
    u = Vector{typeof(u1)}(undef, T)
    u[1] = u1
    y1 = h(p, u[1], 1)
    y = Vector{typeof(y1)}(undef, T)
    y[1] = y1
    for i in 2:T
        u[i] = f(p, u[i - 1], i) .+ g(p, u[i - 1], i) * noise[i]
        y[i] = h(p, u[i], i)
    end
    return y
end

function _solve_zygote(f, g, h, p, u0, noise, tspan, alg::GeneralStateSpace)
    T = tspan[2]
    u1 = f(p, u0, 1) .+ g(p, u0, 1) * noise[1]
    u = Zygote.Buffer(Vector{typeof(u1)}(undef, T))
    u[1] = u1
    y1 = h(p, u[1], 1)
    y = Zygote.Buffer(Vector{typeof(y1)}(undef, T))
    y[1] = y1
    for i in 2:T
        u[i] = f(p, u[i - 1], i) .+ g(p, u[i - 1], i) * noise[i]
        y[i] = h(p, u[i], i)
    end
    copy(u)
    return copy(y)
end

# Linear Time Invariant
function solve(f, g, h, p, u0, noise, tspan, alg::LTIStateSpace)
    T = tspan[2]
    A = f(p)
    B = g(p)
    C = h(p)
    u1 = A * u0 .+ B * noise[1]
    u = Vector{typeof(u1)}(undef, T)
    u[1] = u1
    y1 = C * u[1]
    y = Vector{typeof(y1)}(undef, T)
    y[1] = y1
    for i in 2:T
        u[i] = A * u[i - 1] .+ B * noise[i]
        y[i] = C * u[i]
    end
    return y
end

function _solve_zygote(f, g, h, p, u0, noise, tspan, alg::LTIStateSpace)
    T = tspan[2]
    A = f(p)
    B = g(p)
    C = h(p)
    u1 = A * u0 .+ B * noise[1]
    u = Zygote.Buffer(Vector{typeof(u1)}(undef, T))
    u[1] = u1
    y1 = C * u[1]
    y = Zygote.Buffer(Vector{typeof(y1)}(undef, T))
    y[1] = y1
    for i in 2:T
        u[i] = A * u[i - 1] .+ B * noise[i]
        y[i] = C * u[i]
    end
    copy(u)
    return copy(y)
end

@adjoint function solve(f, g, h, p, u0, noise, tspan, alg)
    return Zygote.pullback(_solve_zygote, f, g, h, p, u0, noise, tspan, alg)
end

##  Example with LTI system.
@adjoint diagm(x::AbstractVector) = diagm(x), dy -> (diag(dy),)
@adjoint function diagm(pr::Pair)
    return diagm(pr), dy -> ((first = nothing, second = diag(dy, first(pr))))
end
# These lines above for diagm can be removed once https://github.com/FluxML/Zygote.jl/issues/570 is fixed
function build_example_LTI(θ, N)
    pad_θ = [θ; zeros(N - length(θ))]
    A = diagm(pad_θ * 0.1)
    B = diagm(pad_θ * 0.01)
    C = diagm(pad_θ * 0.01)
    return ExampleLTI(A, B, C)
end

## Baseline version, nonlinear with time-variation allowed
evolution(p, x, i) = p.A * x
volatility(p, x, i) = p.B
observation(p, x, i) = p.C * x

# LTI versions
evolution(p) = p.A
volatility(p) = p.B
observation(p) = p.C

function test_objective(θ, u0, noise, alg, T, N)
    mod = build_example_LTI(θ, N)
    sol = solve(evolution, volatility, observation, mod, u0, noise, (1, T), alg)  # the model is the `p`
    return mean([sum(vec .^ 2) for vec in sol])
end

# Refresh adjoints!
Zygote.refresh()

## Test starts

# Generate example data
function generate_test_data(n, T, n_θ)
    noise = [randn(n) for _ in 1:T]
    return (θ = randn(n_θ), noise = noise, noise_matrix = hcat(noise...),
            noise_va = VectorOfArray(noise), noise_dea = DiffEqArray(noise, 1:T),  # even slower outside of the primal
            u0 = zeros(n), n = n, T = T, n_θ = n_θ)
end

function benchmark_variations(dat)
    @unpack n, T, θ, n_θ, noise, noise_matrix, noise_va, noise_dea, u0 = dat

    println("Benchmarking n = $n, T = $T\n")
    println("Verify primal speed are all the same")
    alg_general = GeneralStateSpace()
    alg_LTI = LTIStateSpace()

    @btime test_objective($θ, $u0, $noise, $alg_general, T, $n)
    @btime test_objective($θ, $u0, $noise_va, $alg_general, T, $n)
    @btime test_objective($θ, $u0, $noise, $alg_LTI, T, $n)

    println("Gradient | array of arrays | fully nonlinear")
    @btime gradient((θ, u0, noise) -> test_objective(θ, u0, noise, $alg_general, T, $n), $θ,
                    $u0, $noise)
    println("Gradient | VectorOfArrays | fully nonlinear")
    @btime gradient((θ, u0, noise) -> test_objective(θ, u0, noise, $alg_general, T, $n), $θ,
                    $u0, $noise_va)
    println("Gradient | array of arrays | linear-only")
    @btime gradient((θ, u0, noise) -> test_objective(θ, u0, noise, $alg_LTI, T, $n), $θ, $u0,
                    $noise)
    println("Gradient of only θ | array of arrays | fully nonlinear")
    @btime gradient(θ -> test_objective(θ, $u0, $noise, $alg_general, T, $n), $θ)
    println("ForwardDiff Gradient of only θ | array of arrays | fully nonlinear")
    @btime Zygote.forward_jacobian(θ -> test_objective(θ, $u0, $noise, $alg_general, T, $n), $θ)
    return nothing
end

# Run example cases
T = 200
n_θ_1 = 5
n_1 = 5
n_θ_2 = 25
n_2 = 50

dat_1 = generate_test_data(n_1, T, n_θ_1)
dat_2 = generate_test_data(n_2, T, n_θ_2)

benchmark_variations(dat_1)
benchmark_variations(dat_2)

#=
# Check pullbacks
dat = generate_test_data(50, 200, 25)
tmp_mod = build_example_LTI(dat.θ, dat.n)
alg_general = GeneralStateSpace()
alg_LTI = LTIStateSpace()
sol_general, pb_general = Zygote.pullback(solve, evolution, volatility, observation, tmp_mod, dat.u0, dat.noise, (1, 200), alg_general)
sol_LTI, pb_LTI = Zygote.pullback(solve, evolution, volatility, observation, tmp_mod, dat.u0, dat.noise, (1, 200), alg_LTI)
@btime pb_general($sol_general)
@btime pb_LTI($sol_LTI)
=#

## Older variations.  All weakly dominated or intermediate steps in eliminating variations

## Calling the evolution, observation, etc.
## Pass in matrices only
# evolution_mat(A, x, i) = A * x
# observation_mat(C, x, i) = C * x
# volatility_mat(B, x, i) = B

# function solve_matrix_functions(p, u0, noise)
#     T = length(noise)
#     A = p.A
#     B = p.B
#     C = p.C
#     u1 = evolution_mat(A, u0, 1) .+ volatility_mat(B, u0, 1) * noise[1]
#     u = Vector{typeof(u1)}(undef, T)
#     u[1] = u1
#     y1 = observation_mat(C, u[1], 1)
#     y = Vector{typeof(y1)}(undef, T)
#     y[1] = y1
#     for i in 2:T
#         u[i] = evolution_mat(A, u[i - 1], i) .+ volatility_mat(B, u[i - 1], i) * noise[i]
#         y[i] = observation_mat(C, u[i], i)
#     end
#     return y
# end

# function _solve_matrix_functions_zygote(p, u0, noise)
#     T = length(noise)
#     A = p.A
#     B = p.B
#     C = p.C
#     u1 = evolution_mat(A, u0, 1) .+ volatility_mat(B, u0, 1) * noise[1]
#     u = Zygote.Buffer(Vector{typeof(u1)}(undef, T))
#     u[1] = u1
#     y1 = observation_mat(C, u[1], 1)
#     y = Zygote.Buffer(Vector{typeof(y1)}(undef, T))
#     y[1] = y1
#     for i in 2:T
#         u[i] = evolution_mat(A, u[i - 1], i) .+ volatility_mat(B, u[i - 1], i) * noise[i]
#         y[i] = observation_mat(C, u[i], i)
#     end
#     copy(u)
#     return copy(y)
# end

# @adjoint function solve_matrix_functions(p, u0, noise)
#     return Zygote.pullback(_solve_matrix_functions_zygote, p, u0, noise)
# end
# @btime test_objective(solve_matrix_functions, $θ, $u0, $noise)
# println("Gradient | array of arrays | matrix-functions")
# @btime gradient((θ, u0, noise) -> test_objective(solve_matrix_functions, θ, u0, noise), $θ, $u0,
#                 $noise)

## No time-variation, but otherwise fully nonlinear

# evolution_no_t(p, x) = p.A * x
# observation_no_t(p, x) = p.C * x
# volatility_no_t(p, x) = p.B

# function solve_no_t(p, u0, noise)
#     T = length(noise)
#     u1 = evolution_no_t(p, u0) .+ volatility_no_t(p, u0) * noise[1]
#     u = Vector{typeof(u1)}(undef, T)
#     u[1] = u1
#     y1 = observation_no_t(p, u[1])
#     y = Vector{typeof(y1)}(undef, T)
#     y[1] = y1
#     for i in 2:T
#         u[i] = evolution_no_t(p, u[i - 1]) .+ volatility_no_t(p, u[i - 1]) * noise[i]
#         y[i] = observation_no_t(p, u[i])
#     end
#     return y
# end

# function _solve_no_t_zygote(p, u0, noise)
#     T = length(noise)
#     u1 = evolution_no_t(p, u0) .+ volatility_no_t(p, u0) * noise[1]
#     u = Zygote.Buffer(Vector{typeof(u1)}(undef, T))
#     u[1] = u1
#     y1 = observation_no_t(p, u[1])
#     y = Zygote.Buffer(Vector{typeof(y1)}(undef, T))
#     y[1] = y1
#     for i in 2:T
#         u[i] = evolution_no_t(p, u[i - 1]) .+ volatility_no_t(p, u[i - 1]) * noise[i]
#         y[i] = observation_no_t(p, u[i])
#     end
#     copy(u)
#     return copy(y)
# end

# @adjoint solve_no_t(p, u0, noise) = Zygote.pullback(_solve_no_t_zygote, p, u0, noise)

# @btime test_objective(solve_no_t, $θ, $u0, $noise)
# println("Gradient | array of arrays | fully nonlinear, no t")
# @btime gradient((θ, u0, noise) -> test_objective(solve_no_t, θ, u0, noise), $θ, $u0, $noise)
# @btime test_objective(solve_shock_matrix_view_linear_functions, $θ, $u0, $noise_matrix)
# @btime test_objective(solve_shock_matrix_view, $θ, $u0, $noise_va)
# println("Gradient | Matrix | fully nonlinear")
# @btime gradient((θ, u0, noise) -> test_objective(solve_shock_matrix, θ, u0, noise), $θ, $u0,
#                 $noise_matrix)
# println("Gradient | Matrix View | fully nonlinear")
# @btime gradient((θ, u0, noise) -> test_objective(solve_shock_matrix_view, θ, u0, noise), $θ,
#                 $u0, $noise_matrix)
# println("Gradient | Matrix View | linear-only")
# @btime gradient((θ, u0, noise) -> test_objective(solve_shock_matrix_view_linear_functions,
#                                                 θ, u0, noise), $θ, $u0, $noise_matrix)

# ## version where the shocks are a matrix rather than an array of arrays
# function solve_shock_matrix(p, u0, noise::AbstractMatrix)
#     T = size(noise, 2)
#     u1 = evolution(p, u0, 1) .+ volatility(p, u0, 1) * noise[:, 1]
#     u = Vector{typeof(u1)}(undef, T)
#     u[1] = u1
#     y1 = observation(p, u[1], 1)
#     y = Vector{typeof(y1)}(undef, T)
#     y[1] = y1
#     for i in 2:T
#         u[i] = evolution(p, u[i - 1], i) .+ volatility(p, u[i - 1], i) * noise[:, i]
#         y[i] = observation(p, u[i], i)
#     end
#     return y
# end

# function _solve_shock_matrix_zygote(p, u0, noise::AbstractMatrix)
#     T = size(noise, 2)
#     u1 = evolution(p, u0, 1) .+ volatility(p, u0, 1) * noise[:, 1]
#     u = Zygote.Buffer(Vector{typeof(u1)}(undef, T))
#     u[1] = u1
#     y1 = observation(p, u[1], 1)
#     y = Zygote.Buffer(Vector{typeof(y1)}(undef, T))
#     y[1] = y1
#     for i in 2:T
#         u[i] = evolution(p, u[i - 1], i) .+ volatility(p, u[i - 1], i) * noise[:, i]
#         y[i] = observation(p, u[i], i)
#     end
#     copy(u)
#     return copy(y)
# end

# @adjoint function solve_shock_matrix(p, u0, noise)
#     return Zygote.pullback(_solve_shock_matrix_zygote, p, u0, noise)
# end

# function solve_shock_matrix_view(p, u0, noise::AbstractMatrix)
#     T = size(noise, 2)
#     u1 = evolution(p, u0, 1) .+ volatility(p, u0, 1) * view(noise, :, 1)
#     u = Vector{typeof(u1)}(undef, T)
#     u[1] = u1
#     y1 = observation(p, u[1], 1)
#     y = Vector{typeof(y1)}(undef, T)
#     y[1] = y1
#     for i in 2:T
#         u[i] = evolution(p, u[i - 1], i) .+ volatility(p, u[i - 1], i) * view(noise, :, i)
#         y[i] = observation(p, u[i], i)
#     end
#     return y
# end

# function _solve_shock_matrix_view_zygote(p, u0, noise::AbstractMatrix)
#     T = size(noise, 2)
#     u1 = evolution(p, u0, 1) .+ volatility(p, u0, 1) * view(noise, :, 1)
#     u = Zygote.Buffer(Vector{typeof(u1)}(undef, T))
#     u[1] = u1
#     y1 = observation(p, u[1], 1)
#     y = Zygote.Buffer(Vector{typeof(y1)}(undef, T))
#     y[1] = y1
#     for i in 2:T
#         u[i] = evolution(p, u[i - 1], i) .+ volatility(p, u[i - 1], i) * view(noise, :, i)
#         y[i] = observation(p, u[i], i)
#     end
#     copy(u)
#     return copy(y)
# end

# @adjoint function solve_shock_matrix_view(p, u0, noise)
#     return Zygote.pullback(_solve_shock_matrix_view_zygote, p, u0, noise)
# end

# ## Without nested functions and using a matrix.
# function solve_shock_matrix_view_linear_functions(p, u0, noise)
#     T = size(noise, 2)
#     A = p.A
#     B = p.B
#     C = p.C
#     u1 = A * u0 .+ B * view(noise, :, 1)
#     u = Vector{typeof(u1)}(undef, T)
#     u[1] = u1
#     y1 = C * u[1]
#     y = Vector{typeof(y1)}(undef, T)
#     y[1] = y1
#     for i in 2:T
#         u[i] = A * u[i - 1] .+ B * view(noise, :, i)
#         y[i] = C * u[i]
#     end
#     return y
# end

# function _solve_shock_matrix_view_linear_functions_zygote(p, u0, noise)
#     T = size(noise, 2)
#     A = p.A
#     B = p.B
#     C = p.C
#     u1 = A * u0 .+ B * view(noise, :, 1)
#     u = Zygote.Buffer(Vector{typeof(u1)}(undef, T))
#     u[1] = u1
#     y1 = C * u[1]
#     y = Zygote.Buffer(Vector{typeof(y1)}(undef, T))
#     y[1] = y1
#     for i in 2:T
#         u[i] = A * u[i - 1] .+ B * view(noise, :, i)
#         y[i] = C * u[i]
#     end
#     copy(u)
#     return copy(y)
# end

# @adjoint function solve_shock_matrix_view_linear_functions(p, u0, noise)
#     return Zygote.pullback(_solve_shock_matrix_view_linear_functions_zygote, p, u0, noise)
# end

# Without nested functions
# function solve_linear_functions(p, u0, noise)
#     T = length(noise)
#     A = p.A
#     B = p.B
#     C = p.C
#     u1 = A * u0 .+ B * noise[1]
#     u = Vector{typeof(u1)}(undef, T)
#     u[1] = u1
#     y1 = C * u[1]
#     y = Vector{typeof(y1)}(undef, T)
#     y[1] = y1
#     for i in 2:T
#         u[i] = A * u[i - 1] .+ B * noise[i]
#         y[i] = C * u[i]
#     end
#     return y
# end

# function _solve_linear_functions_zygote(p, u0, noise)
#     T = length(noise)
#     A = p.A
#     B = p.B
#     C = p.C
#     u1 = A * u0 .+ B * noise[1]
#     u = Zygote.Buffer(Vector{typeof(u1)}(undef, T))
#     u[1] = u1
#     y1 = C * u[1]
#     y = Zygote.Buffer(Vector{typeof(y1)}(undef, T))
#     y[1] = y1
#     for i in 2:T
#         u[i] = A * u[i - 1] .+ B * noise[i]
#         y[i] = C * u[i]
#     end
#     copy(u)
#     return copy(y)
# end

# @adjoint function solve_linear_functions(p, u0, noise)
#     return Zygote.pullback(_solve_linear_functions_zygote, p, u0, noise)
# end
# @btime test_objective(solve_linear_functions, $θ, $u0, $noise)
# println("Gradient | array of arrays | linear-only")
# @btime gradient((θ, u0, noise) -> test_objective(solve_linear_functions, θ, u0, noise), $θ, $u0,
#                 $noise)
