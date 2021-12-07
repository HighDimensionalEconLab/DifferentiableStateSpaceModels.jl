## General function mapping for the Perturbation Soultions
# ideally, we wouldn't need the specializations for the QTI/etc if these are just as efficient

# First order
dssm_evolution(u, sol::FirstOrderPerturbationSolution, t) = sol.A * u   # f
dssm_volatility(u, sol::FirstOrderPerturbationSolution, t) = sol.B      # g
dssm_observation(u, sol::FirstOrderPerturbationSolution, t) = sol.C * u # h

# # The 2nd order are trickier because of state-space pruning TODO: Move these over to DifferenceEquations.jl
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
