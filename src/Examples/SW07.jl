function SW07()
    # Parameters (with values)
    # ε_w = 10
    # ρ_ga = 0.51
    # ε_p = 10
    # l_bar = 0
    # Π_bar = 0.7
    # B = 0.742
    # μ_w = 0
    # μ_p = 0
    # α = 0.24
    # ψ = 0.2696
    # φ = 6.0144
    # δ = 0.025
    # σ_c = 1.5
    # λ = 0.6361
    # ϕ_p = 1.5
    # ι_w = 0.3243
    # ξ_w = 0.8087
    # ι_p = 0.47
    # ξ_p = 0.6
    # σ_l = 1.9423
    # ϕ_w = 1.5
    # r_π = 1.488
    # r_Δy = 0.2347
    # r_y = 0.0593
    # ρ = 0.8762
    # ρ_a = 0.9977
    # ρ_b = 0.5799
    # ρ_g = 0.9957
    # ρ_i = 0.7165
    # ρ_r = 0
    # ρ_p = 0
    # ρ_w = 0
    # γ_bar = 0.3982
    # gy_ss = 0.18

    ∞ = Inf
    # Parameters
    @variables ε_w, ρ_ga, ε_p, l_bar, Π_bar, B, μ_w, μ_p, α, ψ, φ, δ, σ_c, λ, ϕ_p, ι_w, ξ_w, ι_p, ξ_p, σ_l, ϕ_w, r_π, r_Δy, r_y, ρ, ρ_a, ρ_b, ρ_g, ρ_i, ρ_r, ρ_p, ρ_w, γ_bar, gy_ss # 34 parameters
    # states
    @variables t::Integer, y_f_m(..), y_m(..), k_f(..), k(..), c_f_m(..), c_m(..), i_f_m(..), i_m(..), π_m(..), w_m(..), r_m(..), ε_a(..), b(..), ε_g(..), ε_i(..), ε_r(..), ε_pm(..), η_p_aux(..), ε_wm(..), η_w_aux(..) # 20 variables
    # controls
    @variables k_s_f(..), k_s(..), r_f(..), r(..), r_k_f(..), r_k(..), z_f(..), z(..), w_f(..), w(..), l_f(..), l(..), q_f(..), q(..), y_f(..), y(..), i_f(..), i(..), c_f(..), c(..), π(..), μ_pm(..) # 22 variables

    # Pre-compute parameters
    π_bar = 1 + Π_bar / 100
    γ = 1 + γ_bar / 100
    β = 1 / (1 + B / 100)
    β_bar = β * γ ^ (-σ_c)
    R_star = π_bar / (β * γ ^ (-σ_c))
    R_k_star = 1 / (β * γ^-σ_c) - (1 - δ)
    w_star = (α^α * (1 - α)^(1 - α) / (ϕ_p * R_k_star^α))^(1 / (1 - α))
    k_1 = 1 - (1 - δ) / γ
    i_k = (1 - (1 - δ) / γ) * γ
    l_k = (1 - α) / α * R_k_star / w_star
    k_y = ϕ_p * l_k ^ (α - 1)
    i_y = i_k * k_y
    c_y = 1 - gy_ss - i_k * k_y
    z_y = R_k_star * k_y
    w_hlc = 1 / ϕ_w * (1 - α) / α * R_k_star * k_y / c_y

    # States and controls
    x_sym = [y_f_m, y_m, k_f, k, c_f_m, c_m, i_f_m, i_m, π_m, w_m, r_m, ε_a, b, ε_g, ε_i, ε_r, ε_pm, η_p_aux, ε_wm, η_w_aux]
    y_sym = [k_s_f, k_s, r_f, r, r_k_f, r_k, z_f, z, w_f, w, l_f, l, q_f, q, y_f, y, i_f, i, c_f, c, π, μ_pm]
    # Merge p and p_f
    p = [ε_w, ρ_ga, ε_p, l_bar, Π_bar, B, μ_w, μ_p, α, ψ, φ, δ, σ_c, λ, ϕ_p, ι_w, ξ_w, ι_p, ξ_p, σ_l, ϕ_w, r_π, r_Δy, r_y, ρ, ρ_a, ρ_b, ρ_g, ρ_i, ρ_r, ρ_p, ρ_w, γ_bar, gy_ss]

    # Defining H
    H = [
         α * r_k_f(t) + (1 - α) * w_f(t) - ε_a(t),
         (1 - ψ) / ψ * r_k_f(t) - z_f(t),
         w_f(t) + l_f(t) - k_s_f(t) - r_k_f(t), k_f(t) + z_f(t) - k_s_f(t),
         1 / (1 + β * γ) * (i_f_m(t) + β * γ * i_f(t+1) + 1 / (γ^2 * φ) * q_f(t)) + ε_i(t) - i_f(t),
         -r_f(t) +
         b(t) / ((1 - λ / γ) / (σ_c * (1 + λ / γ))) +
         R_k_star / (R_k_star + 1 - δ) * r_k_f(t+1) +
         (1 - δ) / (R_k_star + 1 - δ) * q_f(t+1) - q_f(t),
         λ / γ / (1 + λ / γ) * c_f_m(t) +
         1 / (1 + λ / γ) * c_f(t+1) +
         (σ_c - 1) * w_hlc / (σ_c * (1 + λ / γ)) * (l_f(t) - l_f(t+1)) -
         (1 - λ / γ) / (σ_c * (1 + λ / γ)) * r_f(t) + b(t) - c_f(t),
         c_y * c_f(t) + i_y * i_f(t) + ε_g(t) + z_y * z_f(t) - y_f(t),
         ϕ_p * (α * k_s_f(t) + (1 - α) * l_f(t) + ε_a(t)) - y_f(t),
         σ_l * l_f(t) + 1 / (1 - λ / γ) * c_f(t) - λ / γ / (1 - λ / γ) * c_f_m(t) - w_f(t),
         (1 - k_1) * k_f(t) + k_1 * i_f(t) + k_1 * γ^2 * φ * ε_i(t) - k_f(t+1),
         # Sticky Part
         α * r_k(t) + (1 - α) * w(t) - ε_a(t) - μ_pm(t),
         (1 - ψ) / ψ * r_k(t) - z(t),
         w(t) + l(t) - k_s(t) - r_k(t), k(t) + z(t) - k_s(t),
         1 / (1 + β * γ) * (i_m(t) + β * γ * i(t+1) + 1 / (γ^2 * φ) * q(t)) + ε_i(t) - i(t),
         -r(t) + π(t+1) +
         b(t) / ((1 - λ / γ) / (σ_c * (1 + λ / γ))) +
         R_k_star / (R_k_star + 1 - δ) * r_k(t+1) +
         (1 - δ) / (R_k_star + 1 - δ) * q(t+1) - q(t),
         λ / γ / (1 + λ / γ) * c_m(t) +
         1 / (1 + λ / γ) * c(t+1) +
         (σ_c - 1) * w_hlc / (σ_c * (1 + λ / γ)) * (l(t) - l(t+1)) -
         (1 - λ / γ) / (σ_c * (1 + λ / γ)) * (r(t) - π(t+1)) + b(t) - c(t),
         c_y * c(t) + i_y * i(t) + ε_g(t) + z_y * z(t) - y(t),
         ϕ_p * (α * k_s(t) + (1 - α) * l(t) + ε_a(t)) - y(t),
         1 / (1 + β * γ * ι_p) * (β * γ * π(t+1) + ι_p * π_m(t) + (1 - ξ_p) * (1 - β * γ * ι_p) / (ξ_p * ((ϕ_p - 1) * ε_p + 1)) * μ_pm(t)) + ε_pm(t) - π(t),
         1 / (1 + β * γ) * w_m(t) + β * γ / (1 + β * γ) * w(t+1) + ι_w / (1 + β * γ) * π_m(t) -
         (1 + β * γ * ι_w) / (1 + β * γ) * π(t) +
         β * γ / (1 + β * γ) * π(t+1) +
         (1 - ξ_w) * (1 - β * γ * ξ_w) / ((1 + β * γ) * ξ_w) / ((ϕ_w - 1) * ε_w + 1) *
         (σ_l * l(t) + 1 / (1 - λ / γ) * c(t) - λ / γ / (1 - λ / γ) * c_m(t) - w(t)) + ε_wm(t) - w(t),
         r_π * (1 - ρ) * π(t) + r_y * (1 - ρ) * (y(t) - y_f(t)) + r_Δy * (y(t) - y_f(t) - y_m(t) + y_f_m(t)) + ρ * r_m(t) + ε_r(t),
         (1 - k_1) * k(t) + k_1 * i(t) + k_1 * γ^2 * φ * ε_i(t) - k(t+1),
         # Shock Transition
         ρ_a * ε_a(t) - ε_a(t+1), ρ_b * b(t) - b(t+1), ρ_g * ε_g(t) - ε_g(t+1), ρ_i * ε_i(t) - ε_i(t+1),
         ρ_r * ε_r(t) - ε_r(t+1), ρ_p * ε_pm(t) + η_p_aux(t+1) - μ_p * η_p_aux(t) - ε_pm(t+1), 0 - η_p_aux(t+1),
         ρ_w * ε_wm(t) + η_w_aux(t+1) - μ_w * η_w_aux(t) - ε_wm(t+1), 0 - η_w_aux(t+1),
         # Variable Mapping
         y_f(t) - y_f_m(t+1), y(t) - y_m(t+1), c_f(t) - c_f_m(t+1), c(t) - c_m(t+1), i_f(t) - i_f_m(t+1), i(t) - i_m(t+1),
         π(t) - π_m(t+1), w(t) - w_m(t+1), r(t) - r_m(t+1)]

    steady_states = [y_f_m(∞) ~ 0, y_m(∞) ~ 0, k_f(∞) ~ 0, k(∞) ~ 0, c_f_m(∞) ~ 0, c_m(∞) ~ 0, i_f_m(∞) ~ 0, i_m(∞) ~ 0, π_m(∞) ~ 0, w_m(∞) ~ 0, r_m(∞) ~ 0, ε_a(∞) ~ 0, b(∞) ~ 0, ε_g(∞) ~ 0, ε_i(∞) ~ 0, ε_r(∞) ~ 0, ε_pm(∞) ~ 0, η_p_aux(∞) ~ 0, ε_wm(∞) ~ 0, η_w_aux(∞) ~ 0, k_s_f(∞) ~ 0, k_s(∞) ~ 0, r_f(∞) ~ 0, r(∞) ~ 0, r_k_f(∞) ~ 0, r_k(∞) ~ 0, z_f(∞) ~ 0, z(∞) ~ 0, w_f(∞) ~ 0, w(∞) ~ 0, l_f(∞) ~ 0, l(∞) ~ 0, q_f(∞) ~ 0, q(∞) ~ 0, y_f(∞) ~ 0, y(∞) ~ 0, i_f(∞) ~ 0, i(∞) ~ 0, c_f(∞) ~ 0, c(∞) ~ 0, π(∞) ~ 0, μ_pm(∞) ~ 0]

    n_ϵ = 2 # TODO: change
    n_x = length(x_sym)
    n_y = length(y_sym)
    n_p = length(p)
    Γ = zeros(Num, n_ϵ, n_ϵ) # TODO: change, also make sure it is not a float64 matrix
    Γ[1, 1] = 1
    Γ[2, 2] = 1
    η = zeros(n_x, n_ϵ) # TODO: change

    return H, (; t, x = x_sym, y = y_sym, p, steady_states, Γ, η, max_order = 1, simplify = false, simplify_Ψ = false, simplify_p = false), "SW07"
end
