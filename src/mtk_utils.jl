# Utilities for working with MTK types for Model construction

match_equation_lhs(var, eq) = (string(var) == string(eq.lhs)) # for ss reordering

# markov variable generation and chaining
make_markov_variable(x) = Num(Variable(Symbol(string(x) * "_p")))
make_ss_variable(x) = Num(Variable(Symbol(string(x) * "_ss")))
function check_markov_name_convention(x)
    x_str = string(x)
    (
        length(x_str) >= 2 &&
        x_str[prevind(x_str, lastindex(x_str)):lastindex(x_str)] == "_p"
    ) && throw(ArgumentError("$x violates the naming convention for a non-Markov variable"))
    return nothing
end
function check_ss_name_convention(x)
    x_str = string(x)
    return (
        length(x_str) >= 3 &&
        x_str[prevind(x_str, lastindex(x_str), 2):lastindex(x_str)] == "_ss"
    ) && throw(
        ArgumentError("$x violates the naming convention for a non-steady-state variable"),
    )
end
function connect_markov_variables(vars)
    check_markov_name_convention.(vars)
    check_ss_name_convention.(vars)
    return [(v, make_markov_variable(v), make_ss_variable(v)) for v in vars]
end
macro make_markov(xs...)
    make_future(x) = Symbol(String(x) * "_p")
    make_ss(x) = Symbol(String(x) * "_ss")
    vars = (make_future.(xs)..., make_ss.(xs)...)
    return esc(ModelingToolkit._parse_vars(:variables, Real, vars))
end

# Takes the vars as a list of variables and sorts the eqs in assignment form: [x ~ a, y ~ b] and sorts by [y, x] etc.
sort_by_variables(eqs, ::Nothing) = eqs
sort_by_variables(::Nothing, vars) = nothing
function sort_by_variables(eqs, vars)
    length(eqs) == length(vars) ||
        throw(ArgumentError("The lengths of the equations and variable list must match"))
    reordered_eqs = Array{Num,1}(undef, 0)
    for var in vars
        matches = filter(eq -> match_equation_lhs(var, eq), eqs)
        length(matches) == 1 ||
            throw(ArgumentError("Each variable must correspond to exactly one equation"))
        push!(reordered_eqs, matches[1].rhs)
    end
    return reordered_eqs
end

## Utilities for Simplifying

substitute_and_simplify(::Nothing, subs) = nothing

function substitute_and_simplify(F::AbstractArray, subs)
    return convert(typeof(F), map(f -> substitute_and_simplify(f, subs), F))
end

function substitute_and_simplify(f, subs)
    return ModelingToolkit.simplify(substitute(f, subs))
end

# Utilities for taking derivatives into the correct structure
recursive_differentiate(::Nothing, ::Nothing) = nothing
recursive_differentiate(::Nothing, x) = nothing
recursive_differentiate(f, ::Nothing) = nothing
function recursive_differentiate(
    f::AbstractVector{Num},
    x::AbstractVector{Num},
)


    F = Num[
        expand_derivatives(
            Differential(ModelingToolkit.value(v))(ModelingToolkit.value(O)),
        ) for O in f, v in x
    ]

    return convert(Matrix{Num}, F)
end
function recursive_differentiate(
    f::AbstractMatrix{Num},
    x::AbstractVector{Num},
)
    # F = [convert(typeof(f), expand_derivatives.(Differential(x_val).(f))) for x_val in x]
    f_value = ModelingToolkit.value(f)
    F = [expand_derivatives.(Differential(ModelingToolkit.value(v)).(f_value)) for v in x]

    return convert(Vector{Matrix{Num}}, F)
end
function recursive_differentiate(
    f::AbstractVector{<:AbstractMatrix{Num}},
    x::AbstractVector{Num},
)
    F = [
        [
            convert(typeof(f_val), expand_derivatives.(Differential(x_val).(f_val))) for
            f_val in f
        ] for x_val in x
    ]
    return convert(Vector{Vector{Matrix{Num}}}, F)
end
function stack_hessians(f::AbstractVector{Num}, x::AbstractVector{Num})
    F = ModelingToolkit.hessian.(f, Ref(x))
    return convert(Vector{Matrix{Num}}, F)
end

# Code to generating function with various
# nothing versions return expressions, which could be evaluated
# note that MTK doesn't support passing "nothing" in for variables, so it converts to [] where appropriate
function build_dssm_function(
    ::Nothing,
    y_p,
    y,
    y_ss,
    x_p,
    x,
    x_ss,
    p,
    p_f)
    return nothing
end
function build_dssm_function(
    ::Nothing,
    y,
    x,
    p,
    p_f)
    return nothing
end
function build_dssm_function(
    ::Nothing,
    p,
    p_f)
    return nothing
end
function build_dssm_function(
    f,
    p,
    p_f)
    return build_function(
        f,
        (isnothing(p) ? [] : p),
        (isnothing(p_f) ? [] : p_f),
        [])[2]
end

function build_dssm_function(
    f,
    y_p,
    y,
    y_ss,
    x_p,
    x,
    x_ss,
    p,
    p_f)
    return build_function(
        f,
        y_p,
        y,
        y_ss,
        x_p,
        x,
        x_ss,
        (isnothing(p) ? [] : p),
        (isnothing(p_f) ? [] : p_f),
        [])[2]
end
function build_dssm_function(
    f,
    y,
    x,
    p,
    p_f)
    return build_function(
        f,
        y,
        x,
        (isnothing(p) ? [] : p),
        (isnothing(p_f) ? [] : p_f),
        [])[2]
end
function build_dssm_function(
    f,
    w,
    p,
    p_f)
    return build_function(
        f,
        w,
        (isnothing(p) ? [] : p),
        (isnothing(p_f) ? [] : p_f),
        []    )[2]
end

# Makes generatlized Generated functions (or nothing)
mk_gg_function(expr) = isnothing(expr) ? nothing : mk_function(Main, expr)
