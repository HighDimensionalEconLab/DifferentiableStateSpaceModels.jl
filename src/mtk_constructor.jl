# Helpers
function default_model_cache_location()
    return joinpath(pkgdir(DifferentiableStateSpaceModels), ".function_cache")
end

function generate_first_order_model(
    H;
    ȳ,
    x̄,
    ȳ_iv,
    x̄_iv,
    Γ,
    Ω,
    η,
    Q,
    p,
    y,
    x,
    p_f,
    verbose,
)

    y_vars = connect_markov_variables(y)
    x_vars = connect_markov_variables(x)
    n_y = length(y_vars)
    n_x = length(x_vars)
    n = n_y + n_x
    n_p = isnothing(p) ? 0 : length(p)
    n_ϵ = size(η, 2)
    n_z = (Q == I) ? n : size(Q, 1)

    # Extract triplets from  y_vars, x_vars for connected variable names
    y = [v[1] for v in y_vars]
    x = [v[1] for v in x_vars]
    y_p = [v[2] for v in y_vars]
    x_p = [v[2] for v in x_vars]
    y_ss = [v[3] for v in y_vars]
    x_ss = [v[3] for v in x_vars]
    y_p_substitutions = [y_p[i] => y[i] for i in eachindex(y)]
    x_p_substitutions = [x_p[i] => x[i] for i in eachindex(x)]
    y_ss_substitutions = [y_ss[i] => y[i] for i in eachindex(y)]
    x_ss_substitutions = [x_ss[i] => x[i] for i in eachindex(x)]
    all_substitutions =
        vcat(y_p_substitutions, x_p_substitutions, y_ss_substitutions, x_ss_substitutions)

    # ensure no reuse of variable
    allunique([p; p_f; y; x; y_p; x_p; y_ss; x_ss]) ||
        throw(ArgumentError("Overlapping variables or parameters"))

    # sort variables in equation assignment form.
    ȳ_iv = sort_by_variables(ȳ_iv, y_ss)
    x̄_iv = sort_by_variables(x̄_iv, x_ss)
    ȳ = sort_by_variables(ȳ, y_ss)
    x̄ = sort_by_variables(x̄, x_ss)

    # steady state requiers differentiation after substitution, and wrt [y; x]
    H̄ = deepcopy(H)
    H̄_sub = substitute_and_simplify(H̄, all_substitutions)

    # Derivatives utilities return nothing if either argument nothing
    H_yp = recursive_differentiate(H, y_p)
    H_y = recursive_differentiate(H, y)
    H_xp = recursive_differentiate(H, x_p)
    H_x = recursive_differentiate(H, x)
    H_p = recursive_differentiate(H, p)
    Γ_p = recursive_differentiate(Γ, p)  # solution object required dense
    Ω_p = recursive_differentiate(Ω, p)
    ȳ_p = recursive_differentiate(ȳ, p)
    x̄_p = recursive_differentiate(x̄, p)
    H̄_w = recursive_differentiate(H̄_sub, [y; x]) # differentiate post-substitution wrt w = [y;x], force a dense one
    H_yp_p = recursive_differentiate(H_yp, p)
    H_xp_p = recursive_differentiate(H_xp, p)
    H_y_p = recursive_differentiate(H_y, p)
    H_x_p = recursive_differentiate(H_x, p)
    Ψ = (n_p == 0) ? nothing : stack_hessians(H, [y_p; y; x_p; x])

    # apply substitutions and simplify if required.
    H_yp_sub = substitute_and_simplify(H_yp, all_substitutions)
    H_xp_sub = substitute_and_simplify(H_xp, all_substitutions)
    H_x_sub = substitute_and_simplify(H_x, all_substitutions)
    H_y_sub = substitute_and_simplify(H_y, all_substitutions)
    H_p_sub = substitute_and_simplify(H_p, all_substitutions)
    H_yp_p_sub = substitute_and_simplify(H_yp_p, all_substitutions)
    H_y_p_sub = substitute_and_simplify(H_y_p, all_substitutions)
    H_xp_p_sub = substitute_and_simplify(H_xp_p, all_substitutions)
    H_x_p_sub = substitute_and_simplify(H_x_p, all_substitutions)
    Ψ_sub = substitute_and_simplify(Ψ, all_substitutions)
    ȳ_p_sub = substitute_and_simplify(ȳ_p, [])
    x̄_p_sub = substitute_and_simplify(x̄_p, [])

    # Generate all functions
    Γ_expr = build_dssm_function(Γ, p, p_f)
    Γ_p_expr = build_dssm_function(Γ_p, p, p_f)
    Ω_expr = build_dssm_function(Ω, p, p_f)
    Ω_p_expr = build_dssm_function(Ω_p, p, p_f)
    H_expr = build_dssm_function(H, y_p, y, y_ss, x_p, x, x_ss, p, p_f)
    H_yp_expr = build_dssm_function(H_yp_sub, y, x, p, p_f)
    H_y_expr = build_dssm_function(H_y_sub, y, x, p, p_f)
    H_xp_expr = build_dssm_function(H_xp_sub, y, x, p, p_f)
    H_x_expr = build_dssm_function(H_x_sub, y, x, p, p_f)
    H_yp_p_expr = build_dssm_function(H_yp_p_sub, y, x, p, p_f)
    H_y_p_expr = build_dssm_function(H_y_p_sub, y, x, p, p_f)
    H_xp_p_expr = build_dssm_function(H_xp_p_sub, y, x, p, p_f)
    H_x_p_expr = build_dssm_function(H_x_p_sub, y, x, p, p_f)
    H_p_expr = build_dssm_function(H_p_sub, y, x, p, p_f)
    Ψ_expr = build_dssm_function(Ψ_sub, y, x, p, p_f)
    H̄_expr = build_dssm_function(H̄_sub, [y; x], p, p_f)
    H̄_w_expr = build_dssm_function(H̄_w, [y; x], p, p_f)
    ȳ_iv_expr = build_dssm_function(ȳ_iv, p, p_f)
    x̄_iv_expr = build_dssm_function(x̄_iv, p, p_f)
    ȳ_expr = build_dssm_function(ȳ, p, p_f)
    x̄_expr = build_dssm_function(x̄, p, p_f)
    ȳ_p_expr = build_dssm_function(ȳ_p_sub, p, p_f)
    x̄_p_expr = build_dssm_function(x̄_p_sub, p, p_f)

    verbose && printstyled("Done Building Model\n", color = :cyan)
    # if the module was included, gets the module name, otherwise returns nothing
    return (;
        n,
        n_y,
        n_x,
        n_p,
        n_ϵ,
        n_z,
        η,
        Q,
        Γ_expr,
        Γ_p_expr,
        Ω_expr,
        Ω_p_expr,
        H_expr,
        H_yp_expr,
        H_y_expr,
        H_xp_expr,
        H_x_expr,
        H_yp_p_expr,
        H_y_p_expr,
        H_xp_p_expr,
        H_x_p_expr,
        H_p_expr,
        Ψ_expr,
        H̄_expr,
        H̄_w_expr,
        ȳ_iv_expr,
        x̄_iv_expr,
        ȳ_expr,
        x̄_expr,
        ȳ_p_expr,
        x̄_p_expr,
    )
end
function save_first_order_module(
    H;
    ȳ = nothing,
    x̄ = nothing,
    ȳ_iv = nothing,
    x̄_iv = nothing,
    Γ,
    Ω = nothing,
    η,
    Q = I,
    p = nothing,
    y,
    x,
    p_f = nothing,
    model_name,
    model_cache_location = default_model_cache_location(),
    overwrite_model_cache = false,
    verbose = false,
)

    model_cache_path = joinpath(model_cache_location, model_name * ".jl")

    # only load cache if the module isn't already loaded in memory
    if (isdefined(Main, Symbol(model_name)) && !overwrite_model_cache)
        verbose && printstyled("Using existing module $model_name\n", color = :cyan)
        return model_cache_path = model_cache_path
    end

    # if path already exists
    if (ispath(model_cache_path) && !overwrite_model_cache)
        # path exists and not overwriting
        verbose &&
            printstyled("Model already generated at $model_cache_path\n", color = :cyan)
    else
        mod = generate_first_order_model(
            H;
            y,
            x,
            ȳ,
            x̄,
            ȳ_iv,
            x̄_iv,
            Γ,
            Ω,
            η,
            Q,
            p,
            p_f,
            verbose,
        )
        @unpack n,
        n_y,
        n_x,
        n_p,
        n_ϵ,
        n_z,
        η,
        Q,
        Γ_expr,
        Γ_p_expr,
        Ω_expr,
        Ω_p_expr,
        H_expr,
        H_yp_expr,
        H_y_expr,
        H_xp_expr,
        H_x_expr,
        H_yp_p_expr,
        H_y_p_expr,
        H_xp_p_expr,
        H_x_p_expr,
        H_p_expr,
        Ψ_expr,
        H̄_expr,
        H̄_w_expr,
        ȳ_iv_expr,
        x̄_iv_expr,
        ȳ_expr,
        x̄_expr,
        ȳ_p_expr,
        x̄_p_expr = mod

        mkpath(model_cache_location)
        open(model_cache_path, "w") do io
            write(io, "module $(model_name)\n")
            write(
                io,
                "using LinearAlgebra, DifferentiableStateSpaceModels, ModelingToolkit\n",
            )
            write(io, "const n_y = $n_y\n")
            write(io, "const n_x = $n_x\n")
            write(io, "const n = $n\n")
            write(io, "const n_p = $n_p\n")
            write(io, "const n_ϵ = $n_ϵ\n")
            write(io, "const n_z = $n_z\n")
            if n_ϵ == 1
                write(io, "const η = reshape($η, $n_x, $n_ϵ)\n")
            else
                write(io, "const η = $η\n")
            end
            write(io, "const Q = $Q\n")
            write(io, "const Γ! = $(Γ_expr)\n")
            write(io, "const Γ_p! = $(Γ_p_expr)\n")
            write(io, "const Ω! = $(Ω_expr)\n")
            write(io, "const Ω_p! = $(Ω_p_expr)\n")
            write(io, "const H! = $(H_expr)\n")
            write(io, "const H_yp! = $(H_yp_expr)\n")
            write(io, "const H_y! = $(H_y_expr)\n")
            write(io, "const H_xp! = $(H_xp_expr)\n")
            write(io, "const H_x! = $(H_x_expr)\n")
            write(io, "const H_yp_p! = $(H_yp_p_expr)\n")
            write(io, "const H_y_p! = $(H_y_p_expr)\n")
            write(io, "const H_xp_p! = $(H_xp_p_expr)\n")
            write(io, "const H_x_p! = $(H_x_p_expr)\n")
            write(io, "const H_p! = $(H_p_expr)\n")
            write(io, "const Ψ! = $(Ψ_expr)\n")
            write(io, "const H̄! = $(H̄_expr)\n")
            write(io, "const H̄_w! = $(H̄_w_expr)\n")
            write(io, "const ȳ_iv! = $(ȳ_iv_expr)\n")
            write(io, "const x̄_iv! = $(x̄_iv_expr)\n")
            write(io, "const ȳ! = $(ȳ_expr)\n")
            write(io, "const x̄! = $(x̄_expr)\n")
            write(io, "const ȳ_p! = $(ȳ_p_expr)\n")
            write(io, "const x̄_p! = $(x̄_p_expr)\n")
            write(io, "const steady_state! = nothing\n")
            return write(io, "end\n") # end module
        end
        verbose && printstyled("Saved $model_name to $model_cache_path\n", color = :cyan)
    end

    # if the module was included, gets the module name, otherwise returns nothing
    return model_cache_path
end


# Probably a cleaner way to avoid the copy/paste, but it is a clear pattern to adapt
function generate_second_order_model(
    H;
    ȳ,
    x̄,
    ȳ_iv,
    x̄_iv,
    Γ,
    Ω,
    η,
    Q,
    p,
    y,
    x,
    p_f,
    verbose,
)

    y_vars = connect_markov_variables(y)
    x_vars = connect_markov_variables(x)
    n_y = length(y_vars)
    n_x = length(x_vars)
    n = n_y + n_x
    n_p = isnothing(p) ? 0 : length(p)
    n_ϵ = size(η, 2)
    n_z = (Q == I) ? n : size(Q, 1)

    # Extract triplets from  y_vars, x_vars for connected variable names
    y = [v[1] for v in y_vars]
    x = [v[1] for v in x_vars]
    y_p = [v[2] for v in y_vars]
    x_p = [v[2] for v in x_vars]
    y_ss = [v[3] for v in y_vars]
    x_ss = [v[3] for v in x_vars]
    y_p_substitutions = [y_p[i] => y[i] for i in eachindex(y)]
    x_p_substitutions = [x_p[i] => x[i] for i in eachindex(x)]
    y_ss_substitutions = [y_ss[i] => y[i] for i in eachindex(y)]
    x_ss_substitutions = [x_ss[i] => x[i] for i in eachindex(x)]
    all_substitutions =
        vcat(y_p_substitutions, x_p_substitutions, y_ss_substitutions, x_ss_substitutions)

    # ensure no reuse of variable
    allunique([p; p_f; y; x; y_p; x_p; y_ss; x_ss]) ||
        throw(ArgumentError("Overlapping variables or parameters"))

    # sort variables in equation assignment form.
    ȳ_iv = sort_by_variables(ȳ_iv, y_ss)
    x̄_iv = sort_by_variables(x̄_iv, x_ss)
    ȳ = sort_by_variables(ȳ, y_ss)
    x̄ = sort_by_variables(x̄, x_ss)

    # steady state requiers differentiation after substitution, and wrt [y; x]
    H̄ = deepcopy(H)
    H̄_sub = substitute_and_simplify(H̄, all_substitutions)

    # Derivatives utilities return nothing if either argument nothing
    H_yp = recursive_differentiate(H, y_p)
    H_y = recursive_differentiate(H, y)
    H_xp = recursive_differentiate(H, x_p)
    H_x = recursive_differentiate(H, x)
    H_p = recursive_differentiate(H, p)
    Γ_p = recursive_differentiate(Γ, p)  # solution object required dense
    Ω_p = recursive_differentiate(Ω, p)
    ȳ_p = recursive_differentiate(ȳ, p)
    x̄_p = recursive_differentiate(x̄, p)
    H̄_w = recursive_differentiate(H̄_sub, [y; x]) # differentiate post-substitution wrt w = [y;x], force a dense one
    H_yp_p = recursive_differentiate(H_yp, p)
    H_xp_p = recursive_differentiate(H_xp, p)
    H_y_p = recursive_differentiate(H_y, p)
    H_x_p = recursive_differentiate(H_x, p)
    Ψ = stack_hessians(H, [y_p; y; x_p; x])
    Ψ_p = (n_p == 0) ? nothing : recursive_differentiate(Ψ, p)
    Ψ_yp = recursive_differentiate(Ψ, y_p)
    Ψ_y = recursive_differentiate(Ψ, y)
    Ψ_xp = recursive_differentiate(Ψ, x_p)
    Ψ_x = recursive_differentiate(Ψ, x)

    # apply substitutions and simplify if required.
    H_yp_sub = substitute_and_simplify(H_yp, all_substitutions)
    H_xp_sub = substitute_and_simplify(H_xp, all_substitutions)
    H_x_sub = substitute_and_simplify(H_x, all_substitutions)
    H_y_sub = substitute_and_simplify(H_y, all_substitutions)
    H_p_sub = substitute_and_simplify(H_p, all_substitutions)
    H_yp_p_sub = substitute_and_simplify(H_yp_p, all_substitutions)
    H_y_p_sub = substitute_and_simplify(H_y_p, all_substitutions)
    H_xp_p_sub = substitute_and_simplify(H_xp_p, all_substitutions)
    H_x_p_sub = substitute_and_simplify(H_x_p, all_substitutions)
    Ψ_sub = substitute_and_simplify(Ψ, all_substitutions)
    Ψ_p_sub = substitute_and_simplify(Ψ_p, all_substitutions)
    Ψ_yp_sub = substitute_and_simplify(Ψ_yp, all_substitutions)
    Ψ_y_sub = substitute_and_simplify(Ψ_y, all_substitutions)
    Ψ_xp_sub = substitute_and_simplify(Ψ_xp, all_substitutions)
    Ψ_x_sub = substitute_and_simplify(Ψ_x, all_substitutions)
    ȳ_p_sub = substitute_and_simplify(ȳ_p, [])
    x̄_p_sub = substitute_and_simplify(x̄_p, [])

    # Generate all functions
    Γ_expr = build_dssm_function(Γ, p, p_f)
    Γ_p_expr = build_dssm_function(Γ_p, p, p_f)
    Ω_expr = build_dssm_function(Ω, p, p_f)
    Ω_p_expr = build_dssm_function(Ω_p, p, p_f)
    H_expr = build_dssm_function(H, y_p, y, y_ss, x_p, x, x_ss, p, p_f)
    H_yp_expr = build_dssm_function(H_yp_sub, y, x, p, p_f)
    H_y_expr = build_dssm_function(H_y_sub, y, x, p, p_f)
    H_xp_expr = build_dssm_function(H_xp_sub, y, x, p, p_f)
    H_x_expr = build_dssm_function(H_x_sub, y, x, p, p_f)
    H_yp_p_expr = build_dssm_function(H_yp_p_sub, y, x, p, p_f)
    H_y_p_expr = build_dssm_function(H_y_p_sub, y, x, p, p_f)
    H_xp_p_expr = build_dssm_function(H_xp_p_sub, y, x, p, p_f)
    H_x_p_expr = build_dssm_function(H_x_p_sub, y, x, p, p_f)
    H_p_expr = build_dssm_function(H_p_sub, y, x, p, p_f)
    Ψ_expr = build_dssm_function(Ψ_sub, y, x, p, p_f)
    Ψ_p_expr = build_dssm_function(Ψ_p_sub, y, x, p, p_f)
    Ψ_yp_expr = build_dssm_function(Ψ_yp_sub, y, x, p, p_f)
    Ψ_y_expr = build_dssm_function(Ψ_y_sub, y, x, p, p_f)
    Ψ_xp_expr = build_dssm_function(Ψ_xp_sub, y, x, p, p_f)
    Ψ_x_expr = build_dssm_function(Ψ_x_sub, y, x, p, p_f)
    H̄_expr = build_dssm_function(H̄_sub, [y; x], p, p_f)
    H̄_w_expr = build_dssm_function(H̄_w, [y; x], p, p_f)
    ȳ_iv_expr = build_dssm_function(ȳ_iv, p, p_f)
    x̄_iv_expr = build_dssm_function(x̄_iv, p, p_f)
    ȳ_expr = build_dssm_function(ȳ, p, p_f)
    x̄_expr = build_dssm_function(x̄, p, p_f)
    ȳ_p_expr = build_dssm_function(ȳ_p_sub, p, p_f)
    x̄_p_expr = build_dssm_function(x̄_p_sub, p, p_f)

    verbose && printstyled("Done Building Model\n", color = :cyan)

    return (;
        n,
        n_y,
        n_x,
        n_p,
        n_ϵ,
        n_z,
        η,
        Q,
        Γ_expr,
        Γ_p_expr,
        Ω_expr,
        Ω_p_expr,
        H_expr,
        H_yp_expr,
        H_y_expr,
        H_xp_expr,
        H_x_expr,
        H_yp_p_expr,
        H_y_p_expr,
        H_xp_p_expr,
        H_x_p_expr,
        H_p_expr,
        Ψ_expr,
        Ψ_p_expr,
        Ψ_yp_expr,
        Ψ_y_expr,
        Ψ_xp_expr,
        Ψ_x_expr,
        H̄_expr,
        H̄_w_expr,
        ȳ_iv_expr,
        x̄_iv_expr,
        ȳ_expr,
        x̄_expr,
        ȳ_p_expr,
        x̄_p_expr,
    )
end
function save_second_order_module(
    H;
    ȳ = nothing,
    x̄ = nothing,
    ȳ_iv = nothing,
    x̄_iv = nothing,
    Γ,
    Ω = nothing,
    η,
    Q = I,
    p = nothing,
    y,
    x,
    p_f = nothing,
    model_name,
    model_cache_location = default_model_cache_location(),
    overwrite_model_cache = false,
    verbose = false,
)

    model_cache_path = joinpath(model_cache_location, model_name * ".jl")

    # only load cache if the module isn't already loaded in memory
    if (isdefined(Main, Symbol(model_name)) && !overwrite_model_cache)
        verbose && printstyled("Using existing module $model_name\n", color = :cyan)
        return model_cache_path = model_cache_path
    end

    # if path already exists
    if (ispath(model_cache_path) && !overwrite_model_cache)
        # path exists and not overwriting
        verbose &&
            printstyled("Model already generated at $model_cache_path\n", color = :cyan)
    else
        mod = generate_second_order_model(
            H;
            y,
            x,
            ȳ,
            x̄,
            ȳ_iv,
            x̄_iv,
            Γ,
            Ω,
            η,
            Q,
            p,
            p_f,
            verbose,
        )
        @unpack n,
        n_y,
        n_x,
        n_p,
        n_ϵ,
        n_z,
        η,
        Q,
        Γ_expr,
        Γ_p_expr,
        Ω_expr,
        Ω_p_expr,
        H_expr,
        H_yp_expr,
        H_y_expr,
        H_xp_expr,
        H_x_expr,
        H_yp_p_expr,
        H_y_p_expr,
        H_xp_p_expr,
        H_x_p_expr,
        H_p_expr,
        Ψ_expr,
        H̄_expr,
        H̄_w_expr,
        ȳ_iv_expr,
        x̄_iv_expr,
        ȳ_expr,
        x̄_expr,
        ȳ_p_expr,
        x̄_p_expr = mod

        @unpack H_yp_p_expr,
        H_y_p_expr,
        H_xp_p_expr,
        H_x_p_expr,
        H_p_expr,
        Ψ_expr,
        Ψ_p_expr,
        Ψ_yp_expr,
        Ψ_y_expr,
        Ψ_xp_expr,
        Ψ_x_expr = mod

        mkpath(model_cache_location)
        open(model_cache_path, "w") do io
            write(io, "module $(model_name)\n")
            write(
                io,
                "using LinearAlgebra, DifferentiableStateSpaceModels, ModelingToolkit\n",
            )
            write(io, "const n_y = $n_y\n")
            write(io, "const n_x = $n_x\n")
            write(io, "const n = $n\n")
            write(io, "const n_p = $n_p\n")
            write(io, "const n_ϵ = $n_ϵ\n")
            write(io, "const n_z = $n_z\n")
            if n_ϵ == 1
                write(io, "const η = reshape($η, $n_x, $n_ϵ)\n")
            else
                write(io, "const η = $η\n")
            end
            write(io, "const Q = $Q\n")
            write(io, "const Γ! = $(Γ_expr)\n")
            write(io, "const Γ_p! = $(Γ_p_expr)\n")
            write(io, "const Ω! = $(Ω_expr)\n")
            write(io, "const Ω_p! = $(Ω_p_expr)\n")
            write(io, "const H! = $(H_expr)\n")
            write(io, "const H_yp! = $(H_yp_expr)\n")
            write(io, "const H_y! = $(H_y_expr)\n")
            write(io, "const H_xp! = $(H_xp_expr)\n")
            write(io, "const H_x! = $(H_x_expr)\n")
            write(io, "const H_yp_p! = $(H_yp_p_expr)\n")
            write(io, "const H_y_p! = $(H_y_p_expr)\n")
            write(io, "const H_xp_p! = $(H_xp_p_expr)\n")
            write(io, "const H_x_p! = $(H_x_p_expr)\n")
            write(io, "const H_p! = $(H_p_expr)\n")
            write(io, "const Ψ! = $(Ψ_expr)\n")
            write(io, "const Ψ_p! = $(Ψ_p_expr)\n")
            write(io, "const Ψ_yp! = $(Ψ_yp_expr)\n")
            write(io, "const Ψ_y! = $(Ψ_y_expr)\n")
            write(io, "const Ψ_xp! = $(Ψ_xp_expr)\n")
            write(io, "const Ψ_x! = $(Ψ_x_expr)\n")
            write(io, "const H̄! = $(H̄_expr)\n")
            write(io, "const H̄_w! = $(H̄_w_expr)\n")
            write(io, "const ȳ_iv! = $(ȳ_iv_expr)\n")
            write(io, "const x̄_iv! = $(x̄_iv_expr)\n")
            write(io, "const ȳ! = $(ȳ_expr)\n")
            write(io, "const x̄! = $(x̄_expr)\n")
            write(io, "const ȳ_p! = $(ȳ_p_expr)\n")
            write(io, "const x̄_p! = $(x̄_p_expr)\n")
            write(io, "const steady_state! = nothing\n")
            return write(io, "end\n") # end module
        end
        verbose && printstyled("Saved $model_name to $model_cache_path\n", color = :cyan)
    end

    # if the module was included, gets the module name, otherwise returns nothing
    return model_cache_path
end
