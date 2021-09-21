function rbc()
    ∞ = Inf
    @variables α, β, ρ, δ, σ
    @variables t::Integer, k(..), z(..), c(..), q(..)
    x = [k, z]
    y = [c, q]
    p = [α, β, ρ, δ, σ]

    H = [
        1 / c(t) - (β / c(t+1)) * (α * exp(z(t+1)) * k(t+1)^(α - 1) + (1 - δ)),
        c(t) + k(t+1) - (1 - δ) * k(t) - q(t),
        q(t) - exp(z(t)) * k(t)^α,
        z(t+1) - ρ * z(t),
    ]

    steady_states = [k(∞) ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
    z(∞) ~ 0,
    c(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
            δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
    q(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)),
    ]

    steady_states_iv = [
    k(∞) ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
    z(∞) ~ 0,
    c(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
            δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
    q(∞) ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)),
    ]

    n_ϵ = 1
    n_x = length(x)
    n_y = length(y)
    n_p = length(p)
    Γ = reshape([σ], n_ϵ, n_ϵ)
    η = reshape([0; -1], n_x, n_ϵ) # η is n_x * n_ϵ matrix

    return H, (; t, x, y, p, steady_states, steady_states_iv, Γ, η), "rbc"
end

function rbc_multiple_shocks()
    H, nt = rbc()
    σ = nt.Γ[1, 1]
    Γ = [σ 0; σ σ]
    η = [0 -1; 0 -1]
    return H, merge(nt, (; Γ, η)), "rbc_multiple_shocks"
end

function rbc_reorder_ss()
    H, nt = rbc()
    steady_states = reverse(nt.steady_states)
    return H, merge(nt, (; steady_states)), "rbc_reorder_ss"
end

function rbc_solve_steady_state()
    H, nt = rbc()
    steady_states = nothing
    return H, merge(nt, (; steady_states)), "rbc_solve_steady_state"
end

function rbc_solve_steady_state_different_iv()
    H, nt = rbc()
    steady_states = nothing

    steady_states_iv = H.steady_states_iv
    steady_states_iv[1].rhs = 0.1
    steady_states_iv[2].rhs = 0.0
    return H, merge(nt, (; steady_states, steady_states_iv)), "rbc_solve_steady_state_different_iv"
end

# Add observation equation with a single source of measurement error
function rbc_observables()
    H, nt = rbc()

    n_x = length(nt.x)
    n_y = length(nt.y)
    n_z = 2 # number of observables
    Q = zeros(n_z, n_x + n_y)
    Q[1, 1] = 1.0
    Q[2, 3] = 1.0

    @variables Ω_1
    Ω = [Ω_1, Ω_1]
    p = [nt.p; Ω_1]
    return H, merge(nt, (; p, Ω, Q)), "rbc_observables"
end

# Add observation equation with two parameters in the source of measurement error
function rbc_observables_separate_variance()
    H, nt = rbc()

    n_x = length(nt.x)
    n_y = length(nt.y)
    n_z = 2 # number of observables
    Q = zeros(n_z, n_x + n_y)
    Q[1, 1] = 1.0
    Q[2, 3] = 1.0

    @variables Ω_1, Ω_2
    Ω = [Ω_1, Ω_2]
    p = [nt.p; Ω_1; Ω_2]
    return H, merge(nt, (; p, Ω, Q)), "rbc_observables_separate_variance"
end