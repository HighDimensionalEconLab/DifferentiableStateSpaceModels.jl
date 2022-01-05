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
struct GeneralLikelihood <: AbstractStateSpaceAlgorithm end
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

function _solve(alg::GeneralLikelihood, noise, observables, f, g, h, D, u0, tspan, p,
                save_everystep, save_posteriors, calculate_logpdf, simulate_observation_noise)

    T = tspan[2]
    @assert tspan[1] == 0
    @assert length(noise) == T

    z0 = h(u0, p, 0)
    u = Vector{typeof(u0)}(undef, T + 1)
    z = Vector{typeof(z0)}(undef, T + 1)
    u[1] = u0
    z[1] = z0
    loglik = 0.0  # remains 0 if no observables
    for t = 2:T+1
        u[t] = f(u[t - 1], p, t - 1) .+ g(u[t - 1], p, t - 1) * noise[t - 1]
        z[t] = h(u[t], p, t)
        loglik += logpdf(D, observables[t - 1] - z[t])  # z_0 doesn't enter likelihood        
    end
    return StateSpaceSolution(nothing, nothing, nothing, nothing, loglik)
end

function ChainRulesCore.rrule(::typeof(_solve), alg::GeneralLikelihood, noise, observables, f, g, h, D, u0, tspan,
                                p::FirstOrderPerturbationSolution,
                                save_everystep, save_posteriors, calculate_logpdf, simulate_observation_noise)
    # Primal computation
    T = tspan[2]

    z0 = h(u0, p, 0)
    g0 = g(u0, p, 0)
    u = Vector{typeof(u0)}(undef, T + 1)
    z = Vector{typeof(z0)}(undef, T + 1)
    g_primal = Vector{typeof(g0)}(undef, T + 1)
    f_pb = Vector{Function}(undef, T + 1)
    g_pb = Vector{Function}(undef, T + 1)
    h_pb = Vector{Function}(undef, T + 1)
    u[1] = u0
    z[1] = z0
    loglik = 0.0  # remains 0 if no observables
    for t in 2:T+1
        # Fetch the pullback functions when iterating the evolution forward
        tmpf, f_pb[t] = Zygote.pullback(f, u[t - 1], p, t - 1)
        g_primal[t], g_pb[t] = Zygote.pullback(g, u[t - 1], p, t - 1)
        u[t] = tmpf + g_primal[t] * noise[t - 1]
        z[t], h_pb[t] = Zygote.pullback(h, u[t], p, t)
        loglik += logpdf(D, observables[t - 1] - z[t])  # z_0 doesn't enter likelihood        
    end
    sol = StateSpaceSolution(nothing, nothing, nothing, nothing, loglik)

    function solve_pb(Δsol)
        Δlogpdf = Δsol.logpdf
        ΔA = similar(p.A)
        ΔB = similar(p.B)
        ΔC = similar(p.C)
        Δnoise = similar(noise)
        Δu = [zero(u0) for _ in 1:T+1]
        fill!(ΔA, 0)
        fill!(ΔB, 0)
        fill!(ΔC, 0)

        for t in (T+1):-1:2
            Δz = -1 * Δlogpdf * Zygote.gradient(logpdf, D, observables[t - 1] - z[t])[2]
            Δh = h_pb[t](Δz)
            Δu[t] += Δh[1]
            Δg = g_pb[t](Δu[t] * noise[t - 1]')
            Δf = f_pb[t](Δu[t])
            Δu[t - 1] = Δf[1] # change that to df + dg
            Δnoise[t - 1] = g_primal[t]' * Δu[t]
            # Now, deal with the coefficients
            ΔA += Δf[2].A
            ΔB += Δg[2].B
            ΔC += Δh[2].C
        end 

        return NoTangent(), NoTangent(), Δnoise, NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(),
        Tangent{typeof(p)}(; A = ΔA, B = ΔB, C = ΔC),
        NoTangent(), NoTangent(), NoTangent(), NoTangent()
    end
    return sol, solve_pb
end

function ChainRulesCore.rrule(::typeof(_solve), alg::GeneralLikelihood, noise, observables, f, g, h, D, u0, tspan,
                                p::SecondOrderPerturbationSolution,
                                save_everystep, save_posteriors, calculate_logpdf, simulate_observation_noise)
    # Primal computation
    T = tspan[2]

    z0 = h(u0, p, 0)
    g0 = g(u0, p, 0)
    u = Vector{typeof(u0)}(undef, T + 1)
    z = Vector{typeof(z0)}(undef, T + 1)
    g_primal = Vector{typeof(g0)}(undef, T + 1)
    f_pb = Vector{Function}(undef, T + 1)
    g_pb = Vector{Function}(undef, T + 1)
    h_pb = Vector{Function}(undef, T + 1)
    u[1] = u0
    z[1] = z0
    loglik = 0.0  # remains 0 if no observables
    for t in 2:T+1
        # Fetch the pullback functions when iterating the evolution forward
        tmpf, f_pb[t] = Zygote.pullback(f, u[t - 1], p, t - 1)
        g_primal[t], g_pb[t] = Zygote.pullback(g, u[t - 1], p, t - 1)
        u[t] = tmpf + g_primal[t] * noise[t - 1]
        z[t], h_pb[t] = Zygote.pullback(h, u[t], p, t)
        loglik += logpdf(D, observables[t - 1] - z[t])  # z_0 doesn't enter likelihood        
    end
    sol = StateSpaceSolution(nothing, nothing, nothing, nothing, loglik)
    
    function solve_pb(Δsol)
        Δlogpdf = Δsol.logpdf
        ΔA_0 = similar(p.A_0)
        ΔA_1 = similar(p.A_1)
        ΔA_2 = similar(p.A_2)
        ΔB = similar(p.B)
        ΔC_0 = similar(p.C_0)
        ΔC_1 = similar(p.C_1)
        ΔC_2 = similar(p.C_2)
        Δnoise = similar(noise)
        Δu = [zero(u0) for _ in 1:T+1]
        fill!(ΔA_0, 0)
        fill!(ΔA_1, 0)
        fill!(ΔA_2, 0)
        fill!(ΔB, 0)
        fill!(ΔC_0, 0)
        fill!(ΔC_1, 0)
        fill!(ΔC_2, 0)

        for t in (T+1):-1:2
            Δz = -1 * Δlogpdf * Zygote.gradient(logpdf, D, observables[t - 1] - z[t])[2]
            Δh = h_pb[t](Δz)
            Δu[t] += Δh[1]
            Δg = g_pb[t](Δu[t] * noise[t - 1]')
            Δf = f_pb[t](Δu[t])
            Δu[t - 1] = Δf[1] # change that to df + dg
            Δnoise[t - 1] = g_primal[t]' * Δu[t]
            # Now, deal with the coefficients
            ΔA_0 += Δf[2].A_0
            ΔA_1 += Δf[2].A_1
            ΔA_2 += Δf[2].A_2
            ΔB += Δg[2].B
            ΔC_0 += Δh[2].C_0
            ΔC_1 += Δh[2].C_1
            ΔC_2 += Δh[2].C_2
        end 

        return NoTangent(), NoTangent(), Δnoise, NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent(),
               Tangent{typeof(p)}(; A_0 = ΔA_0, A_1 = ΔA_1, A_2 = ΔA_2, B = ΔB, C_0 = ΔC_0, C_1 = ΔC_1, C_2 = ΔC_2),
               NoTangent(), NoTangent(), NoTangent(), NoTangent()
    end
    return sol, solve_pb
end

############################################
# LTI Specializations
solve(
    A::AbstractMatrix,
    B::AbstractMatrix,
    C::AbstractMatrix,
    D,
    u0,
    tspan,
    alg = LTI();
    noise = nothing,
    observables = nothing,
) = _solve(alg, noise, observables, A, B, C, D, u0, tspan)

# LTI with observables and without conditioning on noise = Kalman Filter
function _solve(
    alg::LTI,
    noise::Nothing,
    observables,
    A::AbstractMatrix,
    B::AbstractMatrix,
    C::AbstractMatrix,
    D::TuringDiagMvNormal,
    u0,
    tspan,
)
    # hardcoded right now for tspan = (0, T) for T+1 points
    T = tspan[2]
    @assert tspan[1] == 0
    @assert length(observables) == T # i.e. we do not calculate the likelihood of the initial condition

    # Gaussian Prior
    u0_mean = u0.m
    u0_variance = u0.C.U' * u0.C.U
    R = Diagonal(abs2.(D.σ))  # the DistributionsAD doesn't have cov defined the covariance of the MvNormal
    B_prod = B * B'

    # TODO: when saveall = false, etc. don't allocate everything, or at least don't save it
    u = Zygote.Buffer(Vector{Vector{Float64}}(undef, T + 1)) # prior mean
    P = Zygote.Buffer(Vector{Matrix{Float64}}(undef, T + 1)) # prior variance
    z = Zygote.Buffer(Vector{Vector{Float64}}(undef, T + 1)) # mean observation

    u[1] = u0_mean
    P[1] = u0_variance
    z[1] = C * u[1]
    loglik = 0.0

    for i = 2:T+1
        # Kalman iteration
        u[i] = A * u[i-1]
        P[i] = A * P[i-1] * A' + B_prod
        z[i] = C * u[i]

        CP_i = C * P[i]
        V = Symmetric(CP_i * C' + R)
        loglik += logpdf(MvNormal(z[i], V), observables[i-1])
        K = CP_i' / V  # gain
        u[i] += K * (observables[i-1] - z[i])
        P[i] -= K * CP_i
    end
    return StateSpaceSolution(copy(z), copy(u), nothing, copy(P), loglik)
end

# LTI with noise and observables is a joint likelihood
function _solve(
    alg::LTI,
    noise,
    observables,
    A::AbstractMatrix,
    B::AbstractMatrix,
    C::AbstractMatrix,
    D,
    u0,
    tspan,
)
    T = tspan[2]
    @assert tspan[1] == 0
    @assert length(noise) == T

    z0 = C * u0
    u = Zygote.Buffer(Vector{typeof(u0)}(undef, T + 1))
    z = Zygote.Buffer(Vector{typeof(z0)}(undef, T + 1))
    u[1] = u0
    z[1] = z0
    loglik = 0.0  # remains 0 if no observables
    for t = 2:T+1
        u[t] = A * u[t-1] .+ B * noise[t-1]
        z[t] = C * u[t]
        err = observables[t-1] - z[t]
        loglik += logpdf(D, err)  # z_0 doesn't enter likelihood        
    end
    return StateSpaceSolution(copy(z), copy(u), noise, nothing, loglik)
end

# Simulation given the noise and no observables
function _solve(
    alg::LTI,
    noise,
    observables::Nothing,
    A::AbstractMatrix,
    B::AbstractMatrix,
    C::AbstractMatrix,
    D,
    u0,
    tspan,
)
    T = tspan[2]
    @assert tspan[1] == 0
    @assert length(noise) == T

    z0 = C * u0
    u = Zygote.Buffer(Vector{typeof(u0)}(undef, T + 1))
    z = Zygote.Buffer(Vector{typeof(z0)}(undef, T + 1))
    u[1] = u0
    z[1] = z0
    for t = 2:T+1
        u[t] = A * u[t-1] .+ B * noise[t-1]
        z[t] = C * u[t]
    end
    return StateSpaceSolution(copy(z), copy(u), noise, nothing, nothing)
end

# LTILikelihood Only calculate the likelihood
# Hopefully unnecssary with saving options/etc. but putting in to ensure no overhead
function _solve(
    alg::LTILikelihood,
    noise::Nothing,
    observables,
    A::AbstractMatrix,
    B::AbstractMatrix,
    C::AbstractMatrix,
    D::TuringDiagMvNormal,
    u0,
    tspan,
)
    # hardcoded right now for tspan = (0, T) for T+1 points
    T = tspan[2]
    @assert tspan[1] == 0
    @assert length(observables) == T # i.e. we do not calculate the likelihood of the initial condition

    # Gaussian Prior
    u0_mean = u0.m
    u0_variance = u0.C.U' * u0.C.U
    R = Diagonal(abs2.(D.σ))  # the DistributionsAD doesn't have cov defined the covariance of the MvNormal
    B_prod = B * B'

    # TODO: when saveall = false, etc. don't allocate everything, or at least don't save it
    u = Zygote.Buffer(Vector{Vector{Float64}}(undef, T + 1)) # prior mean
    P = Zygote.Buffer(Vector{Matrix{Float64}}(undef, T + 1)) # prior variance
    z = Zygote.Buffer(Vector{Vector{Float64}}(undef, T + 1)) # mean observation

    u[1] = u0_mean
    P[1] = u0_variance
    z[1] = C * u[1]
    loglik = 0.0

    for i = 2:T+1
        # Kalman iteration
        u[i] = A * u[i-1]
        P[i] = A * P[i-1] * A' + B_prod
        z[i] = C * u[i]

        CP_i = C * P[i]
        V = Symmetric(CP_i * C' + R)
        loglik += logpdf(MvNormal(z[i], V), observables[i-1])
        K = CP_i' / V
        # K = cholbackslash(V,CP_i)   # gain
        u[i] += K * (observables[i-1] - z[i])
        P[i] -= K * CP_i
    end
    return StateSpaceSolution(nothing, nothing, nothing, nothing, loglik)
end

function _solve(
    alg::LTILikelihood,
    noise,
    observables,
    A::AbstractMatrix,
    B::AbstractMatrix,
    C::AbstractMatrix,
    D,
    u0,
    tspan,
)
    T = tspan[2]
    @assert tspan[1] == 0
    @assert length(noise) == T

    z0 = C * u0
    u = Zygote.Buffer(Vector{typeof(u0)}(undef, T + 1))
    z = Zygote.Buffer(Vector{typeof(z0)}(undef, T + 1))
    u[1] = u0
    z[1] = z0
    loglik = 0.0  # remains 0 if no observables
    for t = 2:T+1
        u[t] = A * u[t-1] .+ B * noise[t-1]
        z[t] = C * u[t]
        err = observables[t-1] - z[t]
        loglik += logpdf(D, err)  # z_0 doesn't enter likelihood        
    end
    return StateSpaceSolution(nothing, nothing, nothing, nothing, loglik)
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
    solve(sol.A_0, sol.A_1, sol.A_2, sol.B, sol.C_0, sol.C_1, sol.C_2, sol.D, u0, tspan, alg; kwargs...)

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