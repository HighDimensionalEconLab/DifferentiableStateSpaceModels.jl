# Development Environment Setup

## Development Environment
1. Install [Julia](https://julialang.org/downloads/)
2. Install [VS Code](https://code.visualstudio.com/) and its latex extension.  Ssee https://github.com/ubcecon/tutorials/blob/master/vscode.md for more details)
3. If on windows, ensure that you have set your git commandlines to be unix-style
```bash
git config --global core.eol lf
git config --global core.autocrlf false
```
4. Clone the repository.  Either:
   - Install github desktop and choose `Code` drop-down on the github webpage, then choose `Open with Github Desktop`
   - Use `> Git Clone` in vscode and choose `https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl` )
   - Start a vscode Julia terminal and then do  `] dev https://github.com/HighDimensionalEconLab/DifferentiableStateSpaceModels.jl` which will download it into `.julia/dev/https://github.com/DifferentiableStateSpaceModels` folder
   
5. Open the folder in vscode, start up a julia terminal and instantiate all of the packages with `] instantiate`.  If you are not using vscode, then manually start Julia with `--project` option in that folder to ensure you have the project file instantiates.
6. Do `] test` in the julia terminal to run the full regression test.
## Code Standards and Design Principles
- Use the unicode math symbol matching the algebra whenever possible.
- Follow https://github.com/jrevels/YASGuide largely.
    - https://github.com/QuantEcon/lecture-source-jl/blob/master/style.md also is useful, but is intended more for "scripting" code rather than package code.
- TBD: may change the code formatter/style after JuliaFormatter makes progress.