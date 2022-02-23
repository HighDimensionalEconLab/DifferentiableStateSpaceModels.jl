
# Development and Benchmarking

## Setup
One time setup:
1. Setup your environment for [VS Code](https://julia.quantecon.org/software_engineering/tools_editors.html), [github](https://julia.quantecon.org/software_engineering/version_control.html) and [unit testing](https://julia.quantecon.org/software_engineering/testing.html).
2. First start up a Julia repl in vscode this project
3. Activate the global environment with `] activate` instead of the project environment
4. Add in global packages for debugging and benchmarking
```julia
] add BenchmarkTools Infiltrator TestEnv PkgBenchmark
```
5. Activate the benchmarking project
```julia
] activate benchmark
```
6. Connect it the current version of the DSSM package,
```julia
] dev .
```
7. Instantiate all benchmarking dependencies,
```julia
] instantiate
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

If you wish to test the integration with a local version of a package, then after each time you do this step you will need to update the manifest accordingly.  For example, if you have downloaded a `DifferenceEquations.jl` in a parallel folder to this one, then you can do
```
] dev ../DifferenceEquations.jl/.
```
And it will use that version until you restart julia.  As `]test` will only use the pinned version, you can replicate the full regression test in that state with `include("test/runtests.jl")`

### Infiltrator and Debugging

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

Always start with the benchmarks activated, i.e. `] activate benchmark`
A few utilities
```julia
using DifferentiableStateSpaceModels, PkgBenchmark
function save_benchmark(results_file = "baseline")
    data = benchmarkpkg(DifferentiableStateSpaceModels; resultfile = joinpath(pkgdir(DifferentiableStateSpaceModels),"benchmark/$results_file.json"))
    export_markdown(joinpath(pkgdir(DifferentiableStateSpaceModels),"benchmark/trial_$results_file.md"), data)
end
function generate_judgement(new_results, old_results = "baseline", judge_file = "judge")
    return export_markdown(joinpath(pkgdir(DifferentiableStateSpaceModels), "benchmark/$judge_file.md"),
                           judge(PkgBenchmark.readresults(joinpath(pkgdir(DifferentiableStateSpaceModels),
                                                                   "benchmark/$new_results.json")),
                                 PkgBenchmark.readresults(joinpath(pkgdir(DifferentiableStateSpaceModels),
                                                                   "benchmark/$old_results.json"))))
end

```
In your terminal
```julia
save_benchmark("test") # default is "baseline"

# Or manually:
# data = benchmarkpkg(DifferentiableStateSpaceModels; resultfile = joinpath(pkgdir(DifferentiableStateSpaceModels),"benchmark/baseline.json"))
# export_markdown(joinpath(pkgdir(DifferentiableStateSpaceModels),"benchmark/trial.md"), data) # can export as markdown
```

To compare against different parameters or after modifications, load the existing baseline and use the `judge` function to compare

```julia
generate_judgement("test") # defaults to generate_judgement("test", "baseline", "judge")
# Or manually
# data = PkgBenchmark.readresults(joinpath(pkgdir(DifferentiableStateSpaceModels),"benchmark/baseline.json"))
# data_2 = benchmarkpkg(DifferentiableStateSpaceModels, BenchmarkConfig(
#                                             env = Dict("JULIA_NUM_THREADS" => 4, "OPENBLAS_NUM_THREADS" => 1),
#                                             juliacmd = `julia -O3`))
# export_markdown(joinpath(pkgdir(DifferentiableStateSpaceModels),"benchmark/judge.md"), judge(data_2, data))
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
