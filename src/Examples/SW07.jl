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

    # Define symbolic ones
    @parameters ε_w ρ_ga ε_p l_bar Π_bar B μ_w μ_p α ψ φ δ σ_c λ ϕ_p ι_w ξ_w ι_p ξ_p σ_l ϕ_w r_π r_Δy r_y ρ ρ_a ρ_b ρ_g ρ_i ρ_r ρ_p ρ_w γ_bar gy_ss # 34 parameters
    # states
    @variables y_f_m y_m k_f k c_f_m c_m i_f_m i_m π_m w_m r_m ε_a b ε_g ε_i ε_r ε_pm η_p_aux ε_wm η_w_aux # 20 variables
    # controls
    @variables k_s_f k_s r_f r r_k_f r_k z_f z w_f w l_f l q_f q y_f y i_f i c_f c π μ_pm # 22 variables
    # Markovize
    @make_markov y_f_m y_m k_f k c_f_m c_m i_f_m i_m π_m w_m r_m ε_a b ε_g ε_i ε_r ε_pm η_p_aux ε_wm η_w_aux k_s_f k_s r_f r r_k_f r_k z_f z w_f w l_f l q_f q y_f y i_f i c_f c π μ_pm # 42 variables

    x_sym = [
        y_f_m,
        y_m,
        k_f,
        k,
        c_f_m,
        c_m,
        i_f_m,
        i_m,
        π_m,
        w_m,
        r_m,
        ε_a,
        b,
        ε_g,
        ε_i,
        ε_r,
        ε_pm,
        η_p_aux,
        ε_wm,
        η_w_aux,
    ]
    y_sym = [
        k_s_f,
        k_s,
        r_f,
        r,
        r_k_f,
        r_k,
        z_f,
        z,
        w_f,
        w,
        l_f,
        l,
        q_f,
        q,
        y_f,
        y,
        i_f,
        i,
        c_f,
        c,
        π,
        μ_pm,
    ]
    p = [
        ε_w,
        ρ_ga,
        ε_p,
        l_bar,
        Π_bar,
        B,
        μ_w,
        μ_p,
        α,
        ψ,
        φ,
        δ,
        σ_c,
        λ,
        ϕ_p,
        ι_w,
        ξ_w,
        ι_p,
        ξ_p,
        σ_l,
        ϕ_w,
        r_π,
        r_Δy,
        r_y,
        ρ,
        ρ_a,
        ρ_b,
        ρ_g,
        ρ_i,
        ρ_r,
        ρ_p,
        ρ_w,
        γ_bar,
        gy_ss,
    ]
    p_f = []

    # Pre-compute variables
    π_bar = 1 + Π_bar / 100
    γ = 1 + γ_bar / 100
    β = 1 / (1 + B / 100)
    β_bar = β * γ^-σ_c
    R_star = π_bar / (β * γ^-σ_c)
    R_k_star = 1 / (β * γ^-σ_c) - (1 - δ)
    w_star = (α^α * (1 - α)^(1 - α) / (ϕ_p * R_k_star^α))^(1 / (1 - α))
    k_1 = 1 - (1 - δ) / γ
    i_k = (1 - (1 - δ) / γ) * γ
    l_k = (1 - α) / α * R_k_star / w_star
    k_y = ϕ_p * l_k^(α - 1)
    i_y = i_k * k_y
    c_y = 1 - gy_ss - i_k * k_y
    z_y = R_k_star * k_y
    w_hlc = 1 / ϕ_w * (1 - α) / α * R_k_star * k_y / c_y

    # Defining H
    H = [
        α * r_k_f + (1 - α) * w_f - ε_a,
        (1 - ψ) / ψ * r_k_f - z_f,
        w_f + l_f - k_s_f - r_k_f,
        k_f + z_f - k_s_f,
        1 / (1 + β * γ) * (i_f_m + β * γ * i_f_p + 1 / (γ^2 * φ) * q_f) + ε_i - i_f,
        -r_f +
        b / ((1 - λ / γ) / (σ_c * (1 + λ / γ))) +
        R_k_star / (R_k_star + 1 - δ) * r_k_f_p +
        (1 - δ) / (R_k_star + 1 - δ) * q_f_p - q_f,
        λ / γ / (1 + λ / γ) * c_f_m +
        1 / (1 + λ / γ) * c_f_p +
        (σ_c - 1) * w_hlc / (σ_c * (1 + λ / γ)) * (l_f - l_f_p) -
        (1 - λ / γ) / (σ_c * (1 + λ / γ)) * r_f + b - c_f,
        c_y * c_f + i_y * i_f + ε_g + z_y * z_f - y_f,
        ϕ_p * (α * k_s_f + (1 - α) * l_f + ε_a) - y_f,
        σ_l * l_f + 1 / (1 - λ / γ) * c_f - λ / γ / (1 - λ / γ) * c_f_m - w_f,
        (1 - k_1) * k_f + k_1 * i_f + k_1 * γ^2 * φ * ε_i - k_f_p,
        # Sticky Part
        α * r_k + (1 - α) * w - ε_a - μ_pm,
        (1 - ψ) / ψ * r_k - z,
        w + l - k_s - r_k,
        k + z - k_s,
        1 / (1 + β * γ) * (i_m + β * γ * i_p + 1 / (γ^2 * φ) * q) + ε_i - i,
        -r +
        π_p +
        b / ((1 - λ / γ) / (σ_c * (1 + λ / γ))) +
        R_k_star / (R_k_star + 1 - δ) * r_k_p +
        (1 - δ) / (R_k_star + 1 - δ) * q_p - q,
        λ / γ / (1 + λ / γ) * c_m +
        1 / (1 + λ / γ) * c_p +
        (σ_c - 1) * w_hlc / (σ_c * (1 + λ / γ)) * (l - l_p) -
        (1 - λ / γ) / (σ_c * (1 + λ / γ)) * (r - π_p) + b - c,
        c_y * c + i_y * i + ε_g + z_y * z - y,
        ϕ_p * (α * k_s + (1 - α) * l + ε_a) - y,
        1 / (1 + β * γ * ι_p) * (
            β * γ * π_p +
            ι_p * π_m +
            (1 - ξ_p) * (1 - β * γ * ι_p) / (ξ_p * ((ϕ_p - 1) * ε_p + 1)) * μ_pm
        ) + ε_pm - π,
        1 / (1 + β * γ) * w_m + β * γ / (1 + β * γ) * w_p + ι_w / (1 + β * γ) * π_m -
        (1 + β * γ * ι_w) / (1 + β * γ) * π +
        β * γ / (1 + β * γ) * π_p +
        (1 - ξ_w) * (1 - β * γ * ξ_w) / ((1 + β * γ) * ξ_w) / ((ϕ_w - 1) * ε_w + 1) *
        (σ_l * l + 1 / (1 - λ / γ) * c - λ / γ / (1 - λ / γ) * c_m - w) +
        ε_wm - w,
        r_π * (1 - ρ) * π +
        r_y * (1 - ρ) * (y - y_f) +
        r_Δy * (y - y_f - y_m + y_f_m) +
        ρ * r_m +
        ε_r,
        (1 - k_1) * k + k_1 * i + k_1 * γ^2 * φ * ε_i - k_p,
        # Shock Transition
        ρ_a * ε_a - ε_a_p,
        ρ_b * b - b_p,
        ρ_g * ε_g - ε_g_p,
        ρ_i * ε_i - ε_i_p,
        ρ_r * ε_r - ε_r_p,
        ρ_p * ε_pm + η_p_aux_p - μ_p * η_p_aux - ε_pm_p,
        0 - η_p_aux_p,
        ρ_w * ε_wm + η_w_aux_p - μ_w * η_w_aux - ε_wm_p,
        0 - η_w_aux_p,
        # Variable Mapping
        y_f - y_f_m_p,
        y - y_m_p,
        c_f - c_f_m_p,
        c - c_m_p,
        i_f - i_f_m_p,
        i - i_m_p,
        π - π_m_p,
        w - w_m_p,
        r - r_m_p,
    ]

    x̄ = [
        y_f_m_ss ~ 0,
        y_m_ss ~ 0,
        k_f_ss ~ 0,
        k_ss ~ 0,
        c_f_m_ss ~ 0,
        c_m_ss ~ 0,
        i_f_m_ss ~ 0,
        i_m_ss ~ 0,
        π_m_ss ~ 0,
        w_m_ss ~ 0,
        r_m_ss ~ 0,
        ε_a_ss ~ 0,
        b_ss ~ 0,
        ε_g_ss ~ 0,
        ε_i_ss ~ 0,
        ε_r_ss ~ 0,
        ε_pm_ss ~ 0,
        η_p_aux_ss ~ 0,
        ε_wm_ss ~ 0,
        η_w_aux_ss ~ 0,
    ]
    ȳ = [
        k_s_f_ss ~ 0,
        k_s_ss ~ 0,
        r_f_ss ~ 0,
        r_ss ~ 0,
        r_k_f_ss ~ 0,
        r_k_ss ~ 0,
        z_f_ss ~ 0,
        z_ss ~ 0,
        w_f_ss ~ 0,
        w_ss ~ 0,
        l_f_ss ~ 0,
        l_ss ~ 0,
        q_f_ss ~ 0,
        q_ss ~ 0,
        y_f_ss ~ 0,
        y_ss ~ 0,
        i_f_ss ~ 0,
        i_ss ~ 0,
        c_f_ss ~ 0,
        c_ss ~ 0,
        π_ss ~ 0,
        μ_pm_ss ~ 0,
    ]

    x̄_iv = deepcopy(x̄)
    ȳ_iv = deepcopy(ȳ)

    n_ϵ = 2 # TODO: change
    n_x = length(x_sym)
    n_y = length(y_sym)
    n_p = length(p)
    n_p = length(p_f)
    Γ = zeros(ModelingToolkit.Constant, n_ϵ, n_ϵ) # TODO: change, also make sure it is not a float64 matrix
    Γ[1, 1] = 1
    Γ[2, 2] = 1
    η = zeros(n_x, n_ϵ) # TODO: change

    return H, (; x = x_sym, y = y_sym, x̄, ȳ, Γ, η, p_f, p, x̄_iv, ȳ_iv), "SW07"
end
