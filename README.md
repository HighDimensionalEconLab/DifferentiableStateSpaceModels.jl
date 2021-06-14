# DifferentiableStateSpaceModels

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jlperla.github.io/DifferentiableStateSpaceModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jlperla.github.io/DifferentiableStateSpaceModels.jl/dev)
[![Build Status](https://github.com/jlperla/DifferentiableStateSpaceModels.jl/workflows/CI/badge.svg)](https://github.com/jlperla/DifferentiableStateSpaceModels.jl/actions)
[![Coverage](https://codecov.io/gh/jlperla/DifferentiableStateSpaceModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlperla/DifferentiableStateSpaceModels.jl)

# Tensorboard
Install tensorboard via pip or conda.  e.g.
```bash
pip install tensorboard
```
You **do not** need to install tensorflow.

Run code such as the `test/rbc_estimation.jl` code that uses a Turing tensorboard logger, then you can access the logs with
```bash
tensorboard --logdir runs
```
from a terminal inside of the project (e.g. in vscode).
 - It will provide you a link to access to watch the logs which is typically http://localhost:6006/#scalars 
 - Within tensorboard, you typically want to choose the "Show data download links" and "Ignore outliers in chart scaling" as both true, which allows you to easily export the figures/etc.
 - If should be possible to download the runs from a different server or cluster by simply copying the `/runs` folder to your local computer and then using `tensorboard --logdir runs` locally.