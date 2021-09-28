# DifferentiableStateSpaceModels

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HighDimensionalEconLab.github.io/DifferentiableStateSpaceModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HighDimensionalEconLab.github.io/DifferentiableStateSpaceModels.jl/dev)
[![Build Status](https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/workflows/CI/badge.svg)](https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/actions)
[![Coverage](https://codecov.io/gh/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl)

## Development and Benchmarking
See [development.md](development.md)

## Solving Perturbations
TBD

## Derivatives of the Perturbation Solvers

## Example Usage for HMC Estimation

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

After installation you should be able to open the [rbc_example.ipynb](rbc_example.ipynb) example.

For the HMC examples and the `rbc_example.ipnb` above, it will save the results to tensorboard.  While examples are running, you can access the logs with
```bash
tensorboard --logdir runs
```
from a terminal inside of the project (e.g. in vscode).
 - It will provide you a link to access to watch the logs which is typically http://localhost:6006/#scalars 
 - Within tensorboard, you typically want to choose the "Show data download links" and "Ignore outliers in chart scaling" as both true, which allows you to easily export the figures/etc.


