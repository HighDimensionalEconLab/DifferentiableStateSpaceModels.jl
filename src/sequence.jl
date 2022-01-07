# Temporary solution type has superset of all results
struct StateSpaceSolution{T1,T2,T3,T4,T5}
    z::T1 # observables, if relevant
    u::T2 # hidden state, or mean of prior if filtering/estimating
    W::T3 # shocks, if not filtering
    P::T4 # prior variances if filtering, would prefer actually having priors instead
    logpdf::T5  #likelihood of observables
end

## Some algorithm types for dispatching
# LTI = Linear Time Invariant
# QTI = Quadratic Time Invariant?  Not sure term, but mean 2nd order quadratic form
abstract type AbstractStateSpaceAlgorithm end
struct GeneralTimeVarying <: AbstractStateSpaceAlgorithm end
struct LTI <: AbstractStateSpaceAlgorithm end
struct QTI <: AbstractStateSpaceAlgorithm end

# Hopefully temporary algorithm types.  Allows us to optimize a likelihood only solution
# then later we can ensure there is no overhead
struct LTILikelihood <: AbstractStateSpaceAlgorithm end
struct QTILikelihood <: AbstractStateSpaceAlgorithm end

#############################################
## Fully general solve
# Calculates the logpdf if observables passed in.
## Can tweak and eventually check lack of overhead
function solve(
    f,
    g,
    u0,
    tspan,
    p = nothing,
    alg = GeneralTimeVarying();
    h = (u, p, t) -> u,
    D = nothing,
    noise = nothing,
    observables = nothing,
    save_everystep = true,
    save_posteriors = true,
    calculate_logpdf = true,
    simulate_observation_noise = false,
)
    return _solve(
        alg,
        noise,
        observables,
        f,
        g,
        h,
        D,
        u0,
        tspan,
        p,
        save_everystep,
        save_posteriors,
        calculate_logpdf,
        simulate_observation_noise,
    )
end

# Likelhood remains 0 if no obserables.  Need to verify no overhead
maybe_logpdf(observables::Nothing, D, z, i) = 0.0
maybe_logpdf(observables, D, z, i) = logpdf(D, observables[i] .- z)

function _solve(
    alg::GeneralTimeVarying,
    noise,
    observables,
    f,
    g,
    h,
    D,
    u0,
    tspan,
    p,
    save_everystep,
    save_posteriors,
    calculate_logpdf,
    simulate_observation_noise,
)
    T = tspan[2]
    @assert tspan[1] == 0
    @assert length(noise) == T

    z0 = h(u0, p, 0)
    u = Zygote.Buffer(Vector{typeof(u0)}(undef, T + 1))
    z = Zygote.Buffer(Vector{typeof(z0)}(undef, T + 1))
    u[1] = u0
    z[1] = z0
    loglik = 0.0  # remains 0 if no observables
    for t = 2:T+1
        u[t] = f(u[t-1], p, t - 1) .+ g(u[t-1], p, t - 1) * noise[t-1]
        z[t] = h(u[t], p, t)
        loglik += maybe_logpdf(observables, D, z[t], t - 1)  # z_0 doesn't enter likelihood        
    end
    return StateSpaceSolution(copy(z), copy(u), noise, nothing, loglik)
end

############################################
# QTI Specializations
# This is doing state-space pruning, so the u0 will be the size of the underlying state rather than the doubling of it used in the general observation equations
solve(
    A_0,
    A_1,
    A_2,
    B,
    C_0,
    C_1,
    C_2,
    D,
    u0,
    tspan,
    alg = QTI();
    noise = nothing,
    observables = nothing,
) = _solve(alg, noise, observables, A_0, A_1, A_2, B, C_0, C_1, C_2, D, u0, tspan)


# With a symmetric A_2 tensor, there should be a specializaiton for the rrule of the quadratic form
#=
function quad(A_2, x)
    n = size(A_2, 1)
    return map( j -> x' * A_2[j, :, :] * x, 1:n)
end
=#
quad(A::AbstractArray{<:Number,3}, x) =
    @matmul c[l] := sum(j) (@matmul [l, j] := sum(k) A[l, j, k] * x[k]) * x[j]

function ChainRulesCore.rrule(::typeof(quad), A::AbstractArray{<:Number,3}, x)
    res = quad(A, x)
    function quad_pb(Δres)
        ΔA = similar(A)
        Δx = zeros(length(x))
        tmp = x * x'
        n = size(A, 1)
        for i in 1:n
            ΔA[i, :, :] .= tmp .* Δres[i]
            Δx += (A[i, :, :] + A[i, :, :]') * x .* Δres[i]
        end
        return NoTangent(), ΔA, Δx
    end
    return res, quad_pb
end

# QTI with noise and observables is a joint likelihood
function _solve(alg::QTI, noise, observables, A_0, A_1, A_2, B, C_0, C_1, C_2, D, u0, tspan)
    T = tspan[2]
    @assert tspan[1] == 0
    @assert length(noise) == T

    z0 = C_0 + C_1 * u0 + quad(C_2, u0)  # the same u0 used for initial observation
    u_f = Zygote.Buffer(Vector{typeof(u0)}(undef, T + 1))
    u = Zygote.Buffer(Vector{typeof(u0)}(undef, T + 1))
    z = Zygote.Buffer(Vector{typeof(z0)}(undef, T + 1))
    u[1] = u0
    u_f[1] = u0
    z[1] = z0
    loglik = 0.0  # remains 0 if no observables
    for t = 2:T+1
        u_f[t] = A_1 * u_f[t-1] .+ B * noise[t-1]
        u[t] = A_0 + A_1 * u[t-1] + quad(A_2, u_f[t-1]) .+ B * noise[t-1]
        z[t] = C_0 + C_1 * u[t] + quad(C_2, u_f[t])
        err = observables[t-1] - z[t]
        loglik += logpdf(D, err)  # z_0 doesn't enter likelihood        
    end
    return StateSpaceSolution(copy(z), copy(u), noise, nothing, loglik)
end


# QTI with noise and observables is a joint likelihood
function _solve(
    alg::QTILikelihood,
    noise,
    observables,
    A_0,
    A_1,
    A_2,
    B,
    C_0,
    C_1,
    C_2,
    D,
    u0,
    tspan,
)
    T = tspan[2]
    @assert tspan[1] == 0
    @assert length(noise) == T

    z0 = C_0 + C_1 * u0 + quad(C_2, u0)  # the same u0 used for initial observation
    u_f = Zygote.Buffer(Vector{typeof(u0)}(undef, T + 1))
    u = Zygote.Buffer(Vector{typeof(u0)}(undef, T + 1))
    z = Zygote.Buffer(Vector{typeof(z0)}(undef, T + 1))
    u[1] = u0
    u_f[1] = u0
    z[1] = z0
    loglik = 0.0  # remains 0 if no observables
    for t = 2:T+1
        u_f[t] = A_1 * u_f[t-1] .+ B * noise[t-1]
        u[t] = A_0 + A_1 * u[t-1] + quad(A_2, u_f[t-1]) .+ B * noise[t-1]
        z[t] = C_0 + C_1 * u[t] + quad(C_2, u_f[t])
        err = observables[t-1] - z[t]
        loglik += logpdf(D, err)  # z_0 doesn't enter likelihood        
    end
    return StateSpaceSolution(copy(z), copy(u), noise, nothing, loglik)
end


# QTI with noise and observables is a joint likelihood
function _solve(
    alg::QTI,
    noise,
    observables::Nothing,
    A_0,
    A_1,
    A_2,
    B,
    C_0,
    C_1,
    C_2,
    D,
    u0,
    tspan,
)
    T = tspan[2]
    @assert tspan[1] == 0
    @assert length(noise) == T

    z0 = C_0 + C_1 * u0 + quad(C_2, u0)  # the same u0 used for initial observation
    u_f = Zygote.Buffer(Vector{typeof(u0)}(undef, T + 1))
    u = Zygote.Buffer(Vector{typeof(u0)}(undef, T + 1))
    z = Zygote.Buffer(Vector{typeof(z0)}(undef, T + 1))
    u[1] = u0
    u_f[1] = u0
    z[1] = z0
    loglik = 0.0  # remains 0 if no observables
    for t = 2:T+1
        u_f[t] = A_1 * u_f[t-1] .+ B * noise[t-1]
        u[t] = A_0 + A_1 * u[t-1] + quad(A_2, u_f[t-1]) .+ B * noise[t-1]
        z[t] = C_0 + C_1 * u[t] + quad(C_2, u_f[t])
    end
    return StateSpaceSolution(copy(z), copy(u), noise, nothing, loglik)
end


## General function mapping for the Perturbation Soultions
# ideally, we wouldn't need the specializations for the QTI/etc if these are just as efficient

# First order
dssm_evolution(x, sol::FirstOrderPerturbationSolution, t) = sol.A * x
dssm_volatility(x, sol::FirstOrderPerturbationSolution, t) = sol.B
dssm_observation(x, sol::FirstOrderPerturbationSolution, t) = sol.C * x

# The 2nd order are trickier because of state-space pruning
function dssm_evolution(x, sol::SecondOrderPerturbationSolution, t)
    # assume x is stacked as [x_f, x]
    # Try a view?
    x_f_new = sol.A_1 * x[1:sol.n_x]
    x_new = sol.A_0 + sol.A_1 * x[(sol.n_x+1):end] + quad(sol.A_2, x[1:sol.n_x])
    return [x_f_new; x_new]
end

dssm_volatility(x, sol::SecondOrderPerturbationSolution, t) = [sol.B; sol.B]

function dssm_observation(x, sol::SecondOrderPerturbationSolution, t)
    # assume x is stacked as [x_f, x]
    # Try a view?
    return sol.C_0 + sol.C_1 * x[(sol.n_x+1):end] + quad(sol.C_2, x[1:sol.n_x])
end

## Mapping to the Pertubation solutions

# Default algorithm is the likelihood-only for now
solve(sol::FirstOrderPerturbationSolution, u0, tspan, alg = LTILikelihood(); kwargs...) =
    solve(sol.A, sol.B, sol.C, sol.D, u0, tspan, alg; kwargs...)

solve(sol::SecondOrderPerturbationSolution, u0, tspan, alg = QTILikelihood(); kwargs...) =
    solve(
        sol.A_0,
        sol.A_1,
        sol.A_2,
        sol.B,
        sol.C_0,
        sol.C_1,
        sol.C_2,
        sol.D,
        u0,
        tspan,
        alg;
        kwargs...,
    )

############################################
# Utility functions
function save_model_results(file_prefix, path, sol::FirstOrderPerturbationSolution, ϵ, z)
    try
        # For first-order ones, save A, B, C, D
        CSV.write(
            joinpath(path, "$(file_prefix)_A.csv"),
            DataFrame(sol.A),
            writeheader = false,
        )
        CSV.write(
            joinpath(path, "$(file_prefix)_B.csv"),
            DataFrame(sol.B),
            writeheader = false,
        )
        CSV.write(
            joinpath(path, "$(file_prefix)_C.csv"),
            DataFrame(sol.C),
            writeheader = false,
        )
        CSV.write(
            joinpath(path, "$(file_prefix)_D.csv"),
            DataFrame(column_name = sol.D.σ),
            writeheader = false,
        ) # Extract the diagonal
        # In addition, save ϵ and z with different file names
        CSV.write(
            joinpath(path, "$(file_prefix)_noise.csv"),
            DataFrame(transpose(hcat(ϵ...))),
            writeheader = false,
        )
        CSV.write(
            joinpath(path, "$(file_prefix)_observables.csv"),
            DataFrame(transpose(hcat(z...))),
            writeheader = false,
        )
        # Save the ergodic distribution for first order
        CSV.write(
            joinpath(path, "$(file_prefix)_ergodic.csv"),
            DataFrame(sol.x_ergodic.C.L * sol.x_ergodic.C.U),
            writeheader = false,
        )
    catch err
        println("Errors occur when saving the CSV files for first-order")
    end
end

function save_model_results(file_prefix, path, sol::SecondOrderPerturbationSolution, ϵ, z)
    try
        # For second-order ones, save A_0, A_1, A_2, B, C_0, C_1, C_2, D
        # Notice that the B and D objects are the same as those in the first-order case
        CSV.write(
            joinpath(path, "$(file_prefix)_A_0.csv"),
            DataFrame(column_name = sol.A_0),
            writeheader = false,
        ) # Vector
        CSV.write(
            joinpath(path, "$(file_prefix)_A_1.csv"),
            DataFrame(sol.A_1),
            writeheader = false,
        )
        CSV.write(
            joinpath(path, "$(file_prefix)_A_2.csv"),
            DataFrame(column_name = vec(sol.A_2)),
            writeheader = false,
        ) # Save as a vector, reshape it when loading
        CSV.write(
            joinpath(path, "$(file_prefix)_B.csv"),
            DataFrame(sol.B),
            writeheader = false,
        )
        CSV.write(
            joinpath(path, "$(file_prefix)_C_0.csv"),
            DataFrame(column_name = sol.C_0),
            writeheader = false,
        ) # Vector
        CSV.write(
            joinpath(path, "$(file_prefix)_C_1.csv"),
            DataFrame(sol.C_1),
            writeheader = false,
        )
        CSV.write(
            joinpath(path, "$(file_prefix)_C_2.csv"),
            DataFrame(column_name = vec(sol.C_2)),
            writeheader = false,
        ) # Save as a vector, reshape it when loading
        CSV.write(
            joinpath(path, "$(file_prefix)_D.csv"),
            DataFrame(column_name = sol.D.σ),
            writeheader = false,
        ) # Extract the diagonal
        # In addition, save ϵ and z with different file names
        CSV.write(
            joinpath(path, "$(file_prefix)_noise.csv"),
            DataFrame(transpose(hcat(ϵ...))),
            writeheader = false,
        )
        CSV.write(
            joinpath(path, "$(file_prefix)_observables.csv"),
            DataFrame(transpose(hcat(z...))),
            writeheader = false,
        )
    catch err
        println("Errors occur when saving the CSV files for second-order")
    end
end