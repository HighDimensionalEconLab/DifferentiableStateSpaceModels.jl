# Utilities for sequential solutions, likelihoods, etc.

function DifferenceEquations.LinearStateSpaceProblem(sol::AbstractFirstOrderPerturbationSolution,
                                                     x0, tspan; observables = nothing,
                                                     noise = nothing)
    return LinearStateSpaceProblem(sol.A, sol.B, x0, tspan; sol.C,
                                   observables_noise = sol.D,
                                   observables, noise, u0_prior = sol.x_ergodic)
end

function DifferenceEquations.QuadraticStateSpaceProblem(sol::AbstractSecondOrderPerturbationSolution,
                                                        x0, tspan; observables = nothing,
                                                        noise = nothing)
    return QuadraticStateSpaceProblem(sol.A_0, sol.A_1, sol.A_2, sol.B, x0, tspan; sol.C_0,
                                      sol.C_1, sol.C_2,
                                      observables_noise = sol.D,
                                      observables, noise)
end