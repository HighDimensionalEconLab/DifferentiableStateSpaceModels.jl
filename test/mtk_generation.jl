using DifferentiableStateSpaceModels, ModelingToolkit, SparseArrays, LinearAlgebra, Parameters, Test
using DifferentiableStateSpaceModels: substitute_and_simplify, recursive_differentiate,
                                       generate_undef_constructor, sort_by_variables,
                                       sparsify_expression, sort_by_variables,
                                       stack_hessians, build_dssm_function

H, nt = Examples.rbc()
functions_type = DenseFunctions()
parallel = ModelingToolkit.SerialForm()
simplify_exp = true
verbose = true
skipzeros = true
@unpack x, y, x̄, ȳ, Γ, η, p, p_f = nt
ȳ_iv = nothing
x̄_iv = nothing
Q = I
Ω = nothing

is_sparse = (typeof(functions_type) <: SparseFunctions)
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
all_substitutions = vcat(y_p_substitutions, x_p_substitutions, y_ss_substitutions,
                         x_ss_substitutions)

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
H̄_sub = substitute_and_simplify(H̄, all_substitutions; simplify = simplify_exp)

# Derivatives utilities return nothing if either argument nothing
H_yp = recursive_differentiate(H, y_p, is_sparse)
H_y = recursive_differentiate(H, y, is_sparse)
H_xp = recursive_differentiate(H, x_p, is_sparse)
H_x = recursive_differentiate(H, x, is_sparse)
H_p = recursive_differentiate(H, p, is_sparse)
Γ_p = recursive_differentiate(Γ, p, false)  # solution object required dense
Ω_p = recursive_differentiate(Ω, p, false)
ȳ_p = recursive_differentiate(ȳ, p, false)
x̄_p = recursive_differentiate(x̄, p, false)
H̄_w = recursive_differentiate(H̄_sub, [y; x], false) # differentiate post-substitution wrt w = [y;x], force a dense one
H_yp_p = recursive_differentiate(H_yp, p, is_sparse)
H_xp_p = recursive_differentiate(H_xp, p, is_sparse)
H_y_p = recursive_differentiate(H_y, p, is_sparse)
H_x_p = recursive_differentiate(H_x, p, is_sparse)
Ψ = (n_p == 0) ? nothing : stack_hessians(H, [y_p; y; x_p; x], is_sparse)


# apply substitutions and simplify if required.
H_yp_sub = substitute_and_simplify(H_yp, all_substitutions; simplify = simplify_exp)
H_xp_sub = substitute_and_simplify(H_xp, all_substitutions; simplify = simplify_exp)
H_x_sub = substitute_and_simplify(H_x, all_substitutions; simplify = simplify_exp)
H_y_sub = substitute_and_simplify(H_y, all_substitutions; simplify = simplify_exp)
H_p_sub = substitute_and_simplify(H_p, all_substitutions; simplify = simplify_exp)
H_yp_p_sub = substitute_and_simplify(H_yp_p, all_substitutions; simplify = simplify_exp)
H_y_p_sub = substitute_and_simplify(H_y_p, all_substitutions; simplify = simplify_exp)
H_xp_p_sub = substitute_and_simplify(H_xp_p, all_substitutions; simplify = simplify_exp)
H_x_p_sub = substitute_and_simplify(H_x_p, all_substitutions; simplify = simplify_exp)
Ψ_sub = substitute_and_simplify(Ψ, all_substitutions; simplify = simplify_exp)
ȳ_p_sub = substitute_and_simplify(ȳ_p, []; simplify = simplify_exp)
x̄_p_sub = substitute_and_simplify(x̄_p, []; simplify = simplify_exp)

# Generate all functions
Γ_expr = build_dssm_function(Γ, p, p_f; parallel, skipzeros)
Γ_p_expr = build_dssm_function(Γ_p, p, p_f; parallel, skipzeros)
Ω_expr = build_dssm_function(Ω, p, p_f; parallel, skipzeros)
Ω_p_expr = build_dssm_function(Ω_p, p, p_f; parallel, skipzeros)
H_expr = build_dssm_function(H, y_p, y, y_ss, x_p, x, x_ss, p, p_f; parallel, skipzeros)
H_yp_expr = build_dssm_function(H_yp_sub, y, x, p, p_f; parallel, skipzeros)
H_y_expr = build_dssm_function(H_y_sub, y, x, p, p_f; parallel, skipzeros)
H_xp_expr = build_dssm_function(H_xp_sub, y, x, p, p_f; parallel, skipzeros)
H_x_expr = build_dssm_function(H_x_sub, y, x, p, p_f; parallel, skipzeros)
H_yp_p_expr = build_dssm_function(H_yp_p_sub, y, x, p, p_f; parallel, skipzeros)
H_y_p_expr = build_dssm_function(H_y_p_sub, y, x, p, p_f; parallel, skipzeros)
H_xp_p_expr = build_dssm_function(H_xp_p_sub, y, x, p, p_f; parallel, skipzeros)
H_x_p_expr = build_dssm_function(H_x_p_sub, y, x, p, p_f; parallel, skipzeros)
H_p_expr = build_dssm_function(H_p_sub, y, x, p, p_f; parallel, skipzeros)
Ψ_expr = build_dssm_function(Ψ_sub, y, x, p, p_f; parallel, skipzeros)
H̄_expr = build_dssm_function(H̄_sub, [y; x], p, p_f; parallel, skipzeros)
H̄_w_expr = build_dssm_function(H̄_w, [y; x], p, p_f; parallel, skipzeros)
ȳ_iv_expr = build_dssm_function(ȳ_iv, p, p_f; parallel, skipzeros)
x̄_iv_expr = build_dssm_function(x̄_iv, p, p_f; parallel, skipzeros)
ȳ_expr = build_dssm_function(ȳ, p, p_f; parallel, skipzeros)
x̄_expr = build_dssm_function(x̄, p, p_f; parallel, skipzeros)
ȳ_p_expr = build_dssm_function(ȳ_p_sub, p, p_f; parallel, skipzeros)
x̄_p_expr = build_dssm_function(x̄_p_sub, p, p_f; parallel, skipzeros)