# DifferentiableStateSpaceModels

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HighDimensionalEconLab.github.io/DifferentiableStateSpaceModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HighDimensionalEconLab.github.io/DifferentiableStateSpaceModels.jl/dev)
[![Build Status](https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/workflows/CI/badge.svg)](https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/actions)
[![Coverage](https://codecov.io/gh/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl)


# Example Usage for HMC Estimation

With Julia 1.6+, install with `] add DifferentiableStateSpaceModels`  or clone and use directly (e.g. see [Develop](development.md))


To easily install all of the dependencies, go to the folder (i.e. ` pkgdir(DifferentiableStateSpaceModels)`) within a julia REPL and instantiate all package dependencies
```julia 
using DifferentiableStateSpaceModels, Pkg
cd(pkgdir(DifferentiableStateSpaceModels))
] activate; instantiate
```
(or call `Pkg.activate("."); Pkg.instantiate()`).

While not required for the package, some of the examples use Tensorboard, which can be installed with pip or conda (you do not need to install the full tensorflow)
```bash
pip install tensorboard
```

After installation you should be able to open the [rbc_example.ipynb](rbc_example.ipynb) example.  If you did not install jupyter separately, then install IJulia with `] add IJulia` in your package manager, and then open up the example in a terminal with
```julia
using DifferentiableStateSpaceModels, IJulia
jupyterlab(dir=pkgdir(DifferentiableStateSpaceModels))
```
Alternatively, if you have Jupyter installed separately then
```julia
using DifferentiableStateSpaceModels
cd(pkgdir(DifferentiableStateSpaceModels))
; jupyter lab
```

For the HMC examples and the `rbc_example.ipnb` above, it will save the results to tensorboard.  While examples are running, you can access the logs with
```bash
tensorboard --logdir runs
```
from a terminal inside of the project (e.g. in vscode).
 - It will provide you a link to access to watch the logs which is typically http://localhost:6006/#scalars 
 - Within tensorboard, you typically want to choose the "Show data download links" and "Ignore outliers in chart scaling" as both true, which allows you to easily export the figures/etc.


# Development

One time setup:
1. First, setup your environment for [VS Code](https://julia.quantecon.org/software_engineering/tools_editors.html), [github](https://julia.quantecon.org/software_engineering/version_control.html) and [unit testing](https://julia.quantecon.org/software_engineering/testing.html).
2. In your global environment, (i.e. start julia without `--project` or use `]activate` to decactivate the current project) add in
   ```julia
   ] add BenchmarkTools Infiltrator TestEnv JuliaFormatter
   ```
3. Until JuliaFormatter is built into vscode, we need to
   - Install the extension, https://marketplace.visualstudio.com/items?itemName=singularitti.vscode-julia-formatter
   - Open up the vscode settings and add
   ```json
    "[julia]": {
        "editor.defaultFormatter": "singularitti.vscode-julia-formatter"
    },
    "juliaFormatter.style": "yas",
    "juliaFormatter.alwaysUseReturn": true,
    ```
   - You should be able to use the `> Format Document` or `> Format Document With...` choices.  Do not use the builtin vscode julia formatter!  If formatting existing code makes major changes, then you likely didn't set the `yas` style there, which is also available in the settings menu.

## Editing and Debugging code

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
