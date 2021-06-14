# Development Environment Setup

## Setting up Development Enviroment
1. Install [Julia](https://julialang.org/downloads/) and [GitHub Desktop](https://desktop.github.com/) - not strictly required but never hurts to have it!
2. Install `vscode` and follow basic instructions in https://github.com/ubcecon/tutorials/blob/master/vscode.md
  - In particular, https://github.com/ubcecon/tutorials/blob/master/vscode.md#julia, making sure to do the code formatter step.
  - and the git settings in https://github.com/ubcecon/tutorials/blob/master/vscode.md#general-packages-and-setup
3. Clone the repo by either:
  - Clicking on the `Code` then `Open in GitHub Desktop`.
  - Alternatively, you can go `] dev https://github.com/jlperla/DifferentiableStateSpaceModels.jl` in a Julia REPL and it will clone it to the `.julia/dev` folder.
  - If you wanted to, you could then drag that folder back into github desktop.
4. Open it in vscode by right-clicking on the folder it installed to, and then opening a vscode project.
5. Open the Julia repl in vscode  (`Ctrl-Shift-P` and then go `Julia REPL` or something to find it.
6. type `] instantiate` to install all of the packages.  Get coffee.
6. In the REPL run `] test` and it should do the full unit test.

## Formatting code
- Assuming that you setup the code formatter correctly in the above instructions, before checking in any code you should go `Ctrl-Shift-P` and type `Formatter`.


## Code Standards and Design Principles
- Use the unicode math symbol matching the algebra whenever possible.
- Follow https://github.com/jrevels/YASGuide largely.
    - https://github.com/QuantEcon/lecture-source-jl/blob/master/style.md also is useful, but is intended more for "scripting" code rather than package code.
