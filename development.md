
# Development and Benchmarking

## Setup
One time setup:
1. First, setup your environment for [VS Code](https://julia.quantecon.org/software_engineering/tools_editors.html), [github](https://julia.quantecon.org/software_engineering/version_control.html) and [unit testing](https://julia.quantecon.org/software_engineering/testing.html).
2. In your global environment, (i.e. start julia without `--project` or use `]activate` to deactivate the current project) add in
   ```julia
   ] add BenchmarkTools Infiltrator TestEnv PkgBenchmark
   ```

## Formatting Code
Within VS Code, format code before committing with `> Format Document` or `> Format Document With...` choices.  This uses the `.JuliaFormatter.toml` file with [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl).

To always format, you can change check  `Editor: Format on Save`  in the VS Code settings

## Editing and Debugging Code

If you open this folder in VS Code, the `Project.toml` at the root is activated rather than the one in the unit tests.
- The `] test` should work without any chances,
- But to step through individual unit tests which may have test-only dependencies, you can use the `TestEnv` package.  To do this, whenever starting the REPL do
```julia
using TestEnv; TestEnv.activate()
```
At that point, you should be able to edit as if the `test/Project.toml` package was activated.  For example, `include("test/runtests.jl")` should be roughly equivalent to `]test`.  

A useful trick for debugging is with `Infiltrator.jl`. Put in a `@exfiltrate`  in the code, (e.g. inside of a DSSM function) and it pushes all local variables into a global associated with the module.

For example, if `call_the_dssm_function_with_exfiltrate` was a function in the DSSM package with `@exfiltrate` in it, then you could do th following in the REPL or a unit test
```julia
call_the_dssm_function_with_exfiltrate()
@show DifferentiableStateSpaceModels.exfiltrated 
```

Do not forget to remove the exfiltrate when committing to the server though, or performance will plummet.

## Benchmarking
This assumes you are running the repository in VS Code (and hence have the project file activated).  If not, then you will need to activate it accordingly (e.g. `--project` when running Julia).

### Running the Full Benchmarks
In your terminal
```julia 
using DifferentiableStateSpaceModels, PkgBenchmark
data = benchmarkpkg(DifferentiableStateSpaceModels; resultfile = joinpath(pkgdir(DifferentiableStateSpaceModels),"benchmark/baseline.json"))
export_markdown(joinpath(pkgdir(DifferentiableStateSpaceModels),"benchmark/trial.md"), data) # can export as markdown
```

If you wanted to run a comparison with different parameters, or after modifications
```julia
data = PkgBenchmark.readresults(joinpath(pkgdir(DifferentiableStateSpaceModels),"benchmark/baseline.json"))
data_2 = benchmarkpkg(DifferentiableStateSpaceModels, BenchmarkConfig(
                                            env = Dict("JULIA_NUM_THREADS" => 8, "OPENBLAS_NUM_THREADS" => 1),
                                            juliacmd = `julia -O3`))
data_judge = judge(data_2, data)  # compare data_2 vs. data baseline
export_markdown(joinpath(pkgdir(DifferentiableStateSpaceModels),"benchmark/judge.md"), data_judge)  
```
Note: Regenerate the tuning file `tune.json` with `benchmarkpkg(DifferentiableStateSpaceModels, retune = true)`

Also, you can comparing full benchmarks between commit hashes with (although benchmarks need to be in both since it will )
```julia
data = judge(DifferentiableStateSpaceModels, "new-commit-hash", "old-commit-hash")
```

### Running Benchmarks During Development

Rather than the whole PkgBenchmark, you can run the individual benchmarks by either first loading them all up
```julia
using DifferentiableStateSpaceModels
include(joinpath(pkgdir(DifferentiableStateSpaceModels), "benchmark/benchmarks.jl"))
```
And then running individual ones

To use:
- To run part of the benchmarks, you can refer to the global `SUITE`.  For example,
```julia
run(SUITE["rbc"]["rbc_observables"]["first_order"], verbose = true)
```
- Or to get specific statistics such as the median (and using postfix)
```julia
SUITE["rbc"]["rbc_observables"]["first_order"]) |> run |> median
```

To compare between changes, save the results and judge the difference (e.g. median).

For example, with a subset of the suite.  Run it and then save the results
```julia
output_path_old = joinpath(pkgdir(DifferentiableStateSpaceModels), "benchmark/rbc_first_order.json")
BenchmarkTools.save(output_path_old, run(SUITE["rbc"]["rbc_observables"]["first_order"], verbose = true))
```
Now you can reload that stored benchmarking later and compare,
```julia
# Make code change and rerun...
results_new = run(SUITE["rbc"]["rbc_observables"]["first_order"], verbose = true)

#Load to compare to the old one
output_path_old = joinpath(pkgdir(DifferentiableStateSpaceModels), "benchmark/rbc_first_order.json")
results_old = BenchmarkTools.load(output_path_old)[1]

judge_results = judge(median(results_new), median(results_old)) # compare the median/etc.
```
