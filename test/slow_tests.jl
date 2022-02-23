# This includes especially slow tests that shouldn't be run automatically on CI/etc.  In
@time include("sgu.jl")
@time include("FVGQ20.jl")