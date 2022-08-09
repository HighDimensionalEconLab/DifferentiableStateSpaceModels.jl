# Utilities for sequential solutions, likelihoods, etc.
function DifferenceEquations.LinearStateSpaceProblem(sol::AbstractFirstOrderPerturbationSolution,
                                                     x0, tspan; observables = nothing,
                                                     noise = nothing,
                                                     u0_prior_mean = zeros(size(sol.A, 1)),
                                                     u0_prior_var = sol.x_ergodic_var,
                                                     observables_noise = sol.D,
                                                     syms = sol.x_symbols)
    return LinearStateSpaceProblem(sol.A, sol.B, x0, tspan; sol.C,
                                   observables_noise,
                                   observables, noise, u0_prior_mean, u0_prior_var,
                                   syms)
end

function DifferenceEquations.QuadraticStateSpaceProblem(sol::AbstractSecondOrderPerturbationSolution,
                                                        x0, tspan; observables = nothing,
                                                        observables_noise = sol.D,
                                                        noise = nothing,
                                                        syms = sol.x_symbols)
    return QuadraticStateSpaceProblem(sol.A_0, sol.A_1, sol.A_2, sol.B, x0, tspan; sol.C_0,
                                      sol.C_1, sol.C_2,
                                      observables_noise,
                                      observables, noise, syms)
end