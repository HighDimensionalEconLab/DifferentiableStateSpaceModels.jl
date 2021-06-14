using BenchmarkTools, Zygote, LinearAlgebra, Parameters, RecursiveArrayTools
using Zygote: @adjoint

function sim(f, g, p, u0, T)
    u1 = f(p, u0)
    u = Zygote.Buffer(Vector{typeof(u1)}(undef, T))
    u[1] = u1
    for i in 2:T
        u[i] = f(p, u[i - 1]) .+ g(p, u[i-1])
    end
    return copy(u)
end

function sim_linear(f, g, p, u0, T)
    A = f(p)
    b = g(p)
    u1 = A * u0 .+ b
    u = Zygote.Buffer(Vector{typeof(u1)}(undef, T))
    u[1] = u1
    for i in 2:T
        u[i] = A * u[i - 1] .+ b
    end
    return copy(u)
end


struct MyType{T,T2}
    A::T
    b::T2
end
@adjoint function MyType(A, b)
    return MyType(A,b), Δ -> (Δ.A,Δ.b)
end
f(p, u) = p.A * u
f(p) = p.A
g(p, u) = p.b
g(p) = p.b


@adjoint diagm(x::AbstractVector) = diagm(x), dy -> (diag(dy),)
@adjoint function diagm(pr::Pair)
    return diagm(pr), dy -> ((first = nothing, second = diag(dy, first(pr))))
end

function make_model(x)
    return MyType(diagm(x), x)
end

function loss_struct(x, u0, T = 200)
    p = make_model(x)
    u = sim(f, g, p, u0, T)
    return sum(sum(u))
end

function loss_struct_linear(x, u0, T = 200)
    p = make_model(x)
    u = sim_linear(f, g, p, u0, T)
    return sum(sum(u))
end
Zygote.refresh()

N = 100
x = rand(N)
u0 = rand(N)
loss_struct(x, u0)
loss_struct_linear(x, u0)


# test primal
@btime loss_struct($x, $u0)
@btime loss_struct_linear($x, $u0)
@btime gradient(loss_struct, $x, $u0)
@btime gradient(loss_struct_linear, $x, $u0)

#=
# Check Zygote pullbacks
p = make_model(x)
u, pb = Zygote.pullback(sim, f, g, p, u0, 200)
u, pb_linear = Zygote.pullback(sim_linear, f, g, p, u0, 200)
@btime pb($u)
@btime pb_linear($u)
=#
