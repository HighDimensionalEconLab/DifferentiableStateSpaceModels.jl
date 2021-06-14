using DifferentiableStateSpaceModels: vech, inv_vech
const UTILITIES = BenchmarkGroup()
Random.seed!(1234)
const A = rand(10, 10)
const Asym = A' * A
const vech_A = vech(Asym)
UTILITIES["vech"] = @benchmarkable vech($Asym)
UTILITIES["inv_vech"] = @benchmarkable inv_vech($vech_A)

# return for the suite
UTILITIES