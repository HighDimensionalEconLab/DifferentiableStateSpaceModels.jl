function rbc()
    @parameters α β ρ δ σ
    @variables k z c q
    @make_markov k z c q
    x = [k, z]
    y = [c, q]
    p = [α, β]  #spread between theta and p by default, but manipulated later.
    p_f = [ρ, δ, σ]

    H = [1 / c - (β / c_p) * (α * exp(z_p) * k_p^(α - 1) + (1 - δ)),
         c + k_p - (1 - δ) * k - q, q - exp(z) * k^α, z_p - ρ * z]

    x̄ = [k_ss ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1)), z_ss ~ 0]
    ȳ = [c_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
                 δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
          q_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1))]

    x̄_iv = [z_ss ~ 0, k_ss ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1))]
    ȳ_iv = [q_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)),
             c_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
                    δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1))]

    n_ϵ = 1
    n_x = length(x)
    n_y = length(y)
    n_p = length(p)
    Γ = reshape([σ], n_ϵ, n_ϵ)
    η = reshape([0; -1], n_x, n_ϵ) # η is n_x * n_ϵ matrix

    return H, (; x, y, x̄, ȳ, Γ, η, p_f, p, x̄_iv, ȳ_iv), "rbc"
end

# Small named variations on the RBC model
function rbc_empty_p_f()
    H, nt = rbc()
    p = [nt.p; nt.p_f]
    p_f = nothing
    return H, merge(nt, (; p, p_f)), "rbc_empty_p"
end

function rbc_empty_p()
    H, nt = rbc()
    p = nothing
    p_f = [nt.p; nt.p_f]
    return H, merge(nt, (; p, p_f)), "rbc_empty_p"
end

function rbc_empty_p_f_multiple_shocks()
    H, nt = rbc()
    p = [nt.p; nt.p_f]
    p_f = nothing
    σ = nt.Γ[1, 1]
    Γ = [σ 0; σ σ]
    η = [0 -1; 0 -1]
    return H, merge(nt, (; p, p_f, Γ, η)), "rbc_empty_p_f_multiple_shocks"
end

function rbc_reorder_ss()
    H, nt = rbc()
    ȳ = reverse(nt.ȳ)
    return H, merge(nt, (; ȳ)), "rbc_reorder_ss"
end

function rbc_solve_steady_state()
    H, nt = rbc()
    ȳ = nothing
    x̄ = nothing
    ȳ_iv = reverse(nt.ȳ_iv) # not sure why this was necessary?  Prevent old name clashes?
    return H, merge(nt, (; ȳ, x̄, ȳ_iv)), "rbc_solve_steady_state"
end

function rbc_solve_steady_state_different_iv()
    H, nt = rbc()
    ȳ = nothing
    x̄ = nothing

    @variables k z c q
    @make_markov k z c q
    x̄_iv = [z_ss ~ 0, k_ss ~ 0.1]
    return H, merge(nt, (; ȳ, x̄, x̄_iv)), "rbc_solve_steady_state_different_iv"
end

function rbc_observables()
    @parameters α β ρ δ σ Ω_1 Ω_2
    @variables k z c q
    @make_markov k z c q
    x = [k, z]
    y = [c, q]
    p = [α, β, ρ, Ω_1]
    p_f = [δ, σ, Ω_2]

    H = [1 / c - (β / c_p) * (α * exp(z_p) * k_p^(α - 1) + (1 - δ)),
         c + k_p - (1 - δ) * k - q, q - exp(z) * k^α, z_p - ρ * z]

    x̄ = [k_ss ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1)), z_ss ~ 0]
    ȳ = [c_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
                 δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
          q_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1))]

    x̄_iv = [z_ss ~ 0, k_ss ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1))]
    ȳ_iv = [q_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)),
             c_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
                    δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1))]

    n_ϵ = 1
    n_z = 2
    n_x = length(x)
    n_y = length(y)
    n_p = length(p)
    Γ = reshape([σ], n_ϵ, n_ϵ)
    η = reshape([0; -1], n_x, n_ϵ) # η is n_x * n_ϵ matrix

    Q = zeros(n_z, n_x + n_y)
    Q[1, 1] = 1.0
    Q[2, 3] = 1.0

    Ω = [Ω_1, Ω_2]

    # Recall the ordering above:
    # p = [α, β, ρ, Ω_1]
    # p_f = [δ, σ, Ω_2]
    select_p_ss_hash = [1 0 0 0; 0 1 0 0]  # only the α, β trigger a new SS calculation.
    select_p_f_ss_hash = [1 0 0]  # only the δ triggers a new SS calculation
    select_p_perturbation_hash = [1 0 0 0; 0 1 0 0; 0 0 1 0]  # the α, β, ρ trigger things
    select_p_f_perturbation_hash = [1 0 0; 0 1 0]  # in reality, the σ only triggers the perturbation if 2nd order

    return H, (; x, y, x̄, ȳ, Γ, η, p_f, p, x̄_iv, ȳ_iv, Ω, Q, select_p_ss_hash, select_p_f_ss_hash, select_p_perturbation_hash, select_p_f_perturbation_hash), "rbc_observables"
end

function rbc_obervables_iv()
    H, nt = rbc_observables()
    ȳ = nothing
    x̄ = nothing

    @variables k z c q
    @make_markov k z c q
    x̄_iv = [z_ss ~ 0, k_ss ~ 0.1]
    return H, merge(nt, (; ȳ, x̄, x̄_iv)), "rbc_observables_iv"
end

function rbc_observables_benchmark()
    @parameters α β ρ δ σ Ω_1
    @variables k z c q
    @make_markov k z c q
    x = [k, z]
    y = [c, q]
    p = [α, β]
    p_f = [ρ, δ, σ, Ω_1]

    H = [1 / c - (β / c_p) * (α * exp(z_p) * k_p^(α - 1) + (1 - δ)),
         c + k_p - (1 - δ) * k - q, q - exp(z) * k^α, z_p - ρ * z]

    x̄ = [k_ss ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1)), z_ss ~ 0]
    ȳ = [c_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
                 δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1)),
          q_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1))]

    x̄_iv = [z_ss ~ 0, k_ss ~ (((1 / β) - 1 + δ) / α)^(1 / (α - 1))]
    ȳ_iv = [q_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)),
             c_ss ~ (((1 / β) - 1 + δ) / α)^(α / (α - 1)) -
                    δ * (((1 / β) - 1 + δ) / α)^(1 / (α - 1))]

    n_ϵ = 1
    n_z = 2
    n_x = length(x)
    n_y = length(y)
    n_p = length(p)
    Γ = reshape([σ], n_ϵ, n_ϵ)
    η = reshape([0; -1], n_x, n_ϵ) # η is n_x * n_ϵ matrix

    Q = zeros(n_z, n_x + n_y)
    Q[1, 1] = 1.0
    Q[2, 3] = 1.0

    Ω = [Ω_1, Ω_1]

    return H, (; x, y, x̄, ȳ, Γ, η, p_f, p, x̄_iv, ȳ_iv, Ω, Q), "rbc_observables_benchmark"
end
