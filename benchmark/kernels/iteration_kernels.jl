using SparseArrays, LinearAlgebra, StaticArrays, BenchmarkTools, SparseArrays
#using CUDA
const use_cuda = false

## First one.  Relatively naive matlab style
# Note that the @view will be slower for smaller matrices!  Should have a small matrix version
function matlab_style(g_x, h_x, η, Q, x_0, ϵ, z)
    n_x = length(x_0)
    n_y = size(g_x, 1)
    T = size(z, 2)
    x = similar(x_0, n_x, T + 1)
    y = similar(ϵ, n_y, T + 1)
    x[:, 1] .= x_0
    y[:, 1] .= g_x * x_0
    @inbounds for t = 1:T
        x[:, t+1] .= h_x * @view(x[:, t]) .+ η * @view(ϵ[:, t+1])
        y[:, t+1] .= g_x * @view(x[:, t])
    end

    # must be separate loop since in Turing code later
    obs_error = similar(z)
    @inbounds for t = 1:T
        obs_error[:, t] .= @view(z[:, t]) .- Q * [y[:, t]; x[:, t]]
    end
    return (; x, y, obs_error = obs_error)
end

# implement for the staticarrays
function oop_vec(g_x, h_x, Q, x_0, η, ϵ_vec, z_vec)
    n_x = length(x_0)
    n_y = size(g_x, 1)
    T = length(z_vec)
    y_0 = g_x * x_0
    w_0 = [y_0; x_0]

    w = Array{typeof(w_0),1}(undef, 0)
    push!(w, w_0)
    x_val = x_0
    @inbounds for t = 1:T
        x_new = h_x * x_val .+ η * ϵ_vec[t]
        push!(w, [g_x * x_val; x_new])
        x_val = x_new
    end

    # must be separate loop since in Turing code later
    obs_error = Array{eltype(z_vec),1}(undef, 0)
    @inbounds for t = 1:T
        obs = z_vec[t] .- Q * w[t]
        push!(obs_error, obs)
    end
    return (w = w, obs_error = obs_error)
end

# does the ηϵ_vec inplace
function combined_w!(ηϵ, z, f_w, f, h_x, g_x, Q, w_0)
    T = size(ηϵ, 2) - 1
    ηϵ[:, 1] .= w_0  # it is the out

    for t = 1:T
        ηϵ[:, t+1] .= f_w * @view(ηϵ[:, t]) .+ @view(ηϵ[:, t+1])
    end

    # must be separate loop since in Turing code later
    for t = 1:T
        z[:, t] .= @view(z[:, t]) .- Q * @view(ηϵ[:, t])
    end
    return nothing
end

function combined_w_mul!(ηϵ, z, f_w, f, h_x, g_x, Q, w_0)
    T = size(ηϵ, 2) - 1
    ηϵ[:, 1] .= w_0  # it is the out

    for t = 1:T
        mul!(@view(ηϵ[:, t+1]), f_w, @view(ηϵ[:, t]), one(eltype(w_0)), one(eltype(w_0)))  # 5 arg mul!
    end

    # must be separate loop since in Turing code later
    for t = 1:T
        mul!(@view(z[:, t]), Q, @view(ηϵ[:, t]), -one(eltype(w_0)), one(eltype(w_0)))
    end
    return nothing
end

function combined_w_mul_vec!(ηϵ_vec, z_vec, f_w, f, h_x, g_x, Q, w_0)
    T = length(ηϵ_vec) - 1
    n_x = size(h_x, 1)
    ηϵ_vec[1] .= w_0  # it is the out

    for t = 1:T
        mul!(ηϵ_vec[t+1], f_w, ηϵ_vec[t], one(eltype(w_0)), one(eltype(w_0)))  # 5 arg mul!
    end

    # must be separate loop since in Turing code later
    for t = 1:T
        mul!(z_vec[t], Q, ηϵ_vec[t], -one(eltype(w_0)), one(eltype(w_0)))
    end
    return nothing
end

function combined_w_mul_h_x!(ηϵ, z, f_w, f, h_x, g_x, Q, w_0)
    T = size(ηϵ, 2) - 1
    n_y = size(g_x, 1)
    ηϵ[:, 1] .= w_0

    # Try the `h_x` separately inplace
    for t = 1:T
        mul!(@view(ηϵ[(n_y+1):end, t+1]), h_x, @view(ηϵ[(n_y+1):end, t]), 1.0, 1.0)  # 5 arg mul!
        mul!(@view(ηϵ[1:n_y, t]), g_x, @view(ηϵ[(n_y+1):end, t]))
    end

    # must be separate loop since in Turing code later
    for t = 1:T
        mul!(@view(z[:, t]), Q, @view(ηϵ[:, t]), -one(eltype(w_0)), one(eltype(w_0)))
    end
    return nothing
end

function create_test_data(;
    n_x,
    n_y,
    n_ϵ,
    n_z,
    T,
    use_cuda,
    η_sparsity = 0.2,
    use_static = false,
)
    n_x + n_y

    # Model definitinos
    h_x = diagm(rand(n_x))  # ensures it is non-explosive.
    g_x = rand(n_y, n_x)
    η = convert(Array{Float64,2}, sprand(Bool, n_x, n_ϵ, η_sparsity))
    Q = rand(n_z, n)

    # Data and simulations
    ϵ = rand(n_ϵ, T + 1)
    ϵ[:, 1] .= 0.0  # ϵ_0 = 0
    x_0 = zeros(n_x)
    z = rand(n_z, T)  # observables

    # Rearrangement in terms of w for some experiments.
    f = [g_x; h_x]
    f_w = [g_x zeros(n_y, n_y); h_x zeros(n_x, n_y)]  # pads to be a w -> w operator so square.  Not obviously better, though!
    η_w = [η; zeros(n_y, n_ϵ)]
    w_0 = [g_x * x_0; x_0]
    ϵ_vec = [ϵ[:, t] for t = 1:T]
    ηϵ = [η * ϵ; zeros(n_y, T + 1)]
    ηϵ_vec = [ηϵ[:, t] for t = 1:(T+1)]
    z_vec = [z[:, t] for t = 1:T]

    # generate cuda versions
    if use_cuda
        @error("commented out")
        # ηϵ_cuda = CuArray{Float32}(ηϵ)
        # z_cuda = CuArray{Float32}(z)
        # f_w_cuda = CuArray{Float32}(f_w)
        # f_cuda = CuArray{Float32}(f)
        # h_x_cuda = CuArray{Float32}(h_x)
        # g_x_cuda = CuArray{Float32}(g_x)
        # Q_cuda = CuArray{Float32}(Q)
        # w_0_cuda = CuArray{Float32}(w_0)
    else
        ηϵ_cuda = nothing
        z_cuda = nothing
        f_w_cuda = nothing
        f_cuda = nothing
        h_x_cuda = nothing
        g_x_cuda = nothing
        Q_cuda = nothing
        w_0_cuda = nothing
    end

    if use_static
        h_x_static = SMatrix{n_x,n_x}(h_x)
        g_x_static = SMatrix{n_y,n_x}(g_x)
        Q_static = SMatrix{n_x,n_x + n_y}(Q)
        x_0_static = SVector{n_x}(x_0)
        ϵ_vec_static = [SVector{n_ϵ}(ϵ_vec[t]) for t = 1:T]
        z_vec_static = [SVector{n_z}(z_vec[t]) for t = 1:T]
        η_static = SMatrix{n_x,n_ϵ}(η)
    else
        h_x_static = nothing
        g_x_static = nothing
        Q_static = nothing
        x_0_static = nothing
        ϵ_vec_static = nothing
        z_vec_static = nothing
        η_static = nothing
    end

    return (;
        h_x,
        g_x,
        η,
        Q,
        ϵ,
        x_0,
        z,
        f,
        f_w,
        η_w,
        w_0,
        ηϵ,
        ηϵ_vec,
        z_vec,
        ϵ_vec,
        ηϵ_cuda,
        z_cuda,
        f_w_cuda,
        f_cuda,
        h_x_cuda,
        g_x_cuda,
        Q_cuda,
        w_0_cuda,
        h_x_static,
        g_x_static,
        Q_static,
        x_0_static,
        ϵ_vec_static,
        z_vec_static,
        η_static,
    )
end

const large = create_test_data(;
    n_x = 40,
    n_y = 30,
    n_ϵ = 20,
    n_z = 20,
    T = 500,
    η_sparsity = 0.2,
    use_cuda = use_cuda,
)
const medium = create_test_data(;
    n_x = 10,
    n_y = 8,
    n_ϵ = 7,
    n_z = 11,
    T = 500,
    η_sparsity = 0.2,
    use_cuda = use_cuda,
)
const small = create_test_data(;
    n_x = 2,
    n_y = 1,
    n_ϵ = 1,
    n_z = 2,
    T = 500,
    η_sparsity = 1.0,
    use_cuda = use_cuda,
    use_static = true,
)

const KERNELS = BenchmarkGroup()
const KERNELS["large"] = BenchmarkGroup()
const KERNELS["medium"] = BenchmarkGroup()
const KERNELS["small"] = BenchmarkGroup()

KERNELS["large"]["matlab_style"] = @benchmarkable matlab_style(
    $large.g_x,
    $large.h_x,
    $large.η,
    $large.Q,
    $large.x_0,
    $large.ϵ,
    $large.z,
)
KERNELS["large"]["combined_w"] = @benchmarkable combined_w!(
    ηϵ,
    z,
    $large.f_w,
    $large.f,
    $large.h_x,
    $large.g_x,
    $large.Q,
    $large.w_0,
) setup = (ηϵ = copy(large.ηϵ);
z = copy(large.z))

KERNELS["large"]["combined_w_mul"] = @benchmarkable combined_w_mul!(
    ηϵ,
    z,
    $large.f_w,
    $large.f,
    $large.h_x,
    $large.g_x,
    $large.Q,
    $large.w_0,
) setup = (ηϵ = copy(large.ηϵ);
z = copy(large.z))
KERNELS["large"]["combined_w_mul_vec"] = @benchmarkable combined_w_mul_vec!(
    ηϵ_vec,
    z_vec,
    $large.f_w,
    $large.f,
    $large.h_x,
    $large.g_x,
    $large.Q,
    $large.w_0,
) setup = (ηϵ_vec = deepcopy(large.ηϵ_vec);
z_vec = deepcopy(large.z_vec))
KERNELS["large"]["combined_w_mul_h_x"] = @benchmarkable combined_w_mul_h_x!(
    ηϵ,
    z,
    $large.f_w,
    $large.f,
    $large.h_x,
    $large.g_x,
    $large.Q,
    $large.w_0,
) setup = (ηϵ = copy(large.ηϵ);
z = copy(large.z))

KERNELS["medium"]["matlab_style"] = @benchmarkable matlab_style(
    $medium.g_x,
    $medium.h_x,
    $medium.η,
    $medium.Q,
    $medium.x_0,
    $medium.ϵ,
    $medium.z,
)
KERNELS["medium"]["combined_w"] = @benchmarkable combined_w!(
    ηϵ,
    z,
    $medium.f_w,
    $medium.f,
    $medium.h_x,
    $medium.g_x,
    $medium.Q,
    $medium.w_0,
) setup = (ηϵ = copy(medium.ηϵ);
z = copy(medium.z))
KERNELS["medium"]["combined_w_mul"] = @benchmarkable combined_w_mul!(
    ηϵ,
    z,
    $medium.f_w,
    $medium.f,
    $medium.h_x,
    $medium.g_x,
    $medium.Q,
    $medium.w_0,
) setup = (ηϵ = copy(medium.ηϵ);
z = copy(medium.z))
KERNELS["medium"]["combined_w_mul_vec"] = @benchmarkable combined_w_mul_vec!(
    ηϵ_vec,
    z_vec,
    $medium.f_w,
    $medium.f,
    $medium.h_x,
    $medium.g_x,
    $medium.Q,
    $medium.w_0,
) setup = (ηϵ_vec = deepcopy(medium.ηϵ_vec);
z_vec = deepcopy(medium.z_vec))
KERNELS["medium"]["combined_w_mul_h_x"] = @benchmarkable combined_w_mul_h_x!(
    ηϵ,
    z,
    $medium.f_w,
    $medium.f,
    $medium.h_x,
    $medium.g_x,
    $medium.Q,
    $medium.w_0,
) setup = (ηϵ = copy(medium.ηϵ);
z = copy(medium.z))
KERNELS["medium"]["oop_vec"] = @benchmarkable oop_vec(
    $medium.g_x,
    $medium.h_x,
    $medium.Q,
    $medium.x_0,
    $medium.η,
    $medium.ϵ_vec,
    $medium.z_vec,
)

KERNELS["small"]["matlab_style"] = @benchmarkable matlab_style(
    $small.g_x,
    $small.h_x,
    $small.η,
    $small.Q,
    $small.x_0,
    $small.ϵ,
    $small.z,
)
KERNELS["small"]["combined_w"] = @benchmarkable combined_w!(
    ηϵ,
    z,
    $small.f_w,
    $small.f,
    $small.h_x,
    $small.g_x,
    $small.Q,
    $small.w_0,
) setup = (ηϵ = copy(small.ηϵ);
z = copy(small.z))
KERNELS["small"]["combined_w_mul"] = @benchmarkable combined_w_mul!(
    ηϵ,
    z,
    $small.f_w,
    $small.f,
    $small.h_x,
    $small.g_x,
    $small.Q,
    $small.w_0,
) setup = (ηϵ = copy(small.ηϵ);
z = copy(small.z))
KERNELS["small"]["combined_w_mul_vec"] = @benchmarkable combined_w_mul_vec!(
    ηϵ_vec,
    z_vec,
    $small.f_w,
    $small.f,
    $small.h_x,
    $small.g_x,
    $small.Q,
    $small.w_0,
) setup = (ηϵ_vec = deepcopy(small.ηϵ_vec);
z_vec = deepcopy(small.z_vec))
KERNELS["small"]["combined_w_mul_h_x"] = @benchmarkable combined_w_mul_h_x!(
    ηϵ,
    z,
    $small.f_w,
    $small.f,
    $small.h_x,
    $small.g_x,
    $small.Q,
    $small.w_0,
) setup = (ηϵ = copy(small.ηϵ);
z = copy(small.z))
KERNELS["small"]["oop_vec"] = @benchmarkable oop_vec(
    $small.g_x,
    $small.h_x,
    $small.Q,
    $small.x_0,
    $small.η,
    $small.ϵ_vec,
    $small.z_vec,
)
KERNELS["small"]["oop_vec_static"] = @benchmarkable oop_vec(
    $small.g_x_static,
    $small.h_x_static,
    $small.Q_static,
    $small.x_0_static,
    $small.η_static,
    $small.ϵ_vec_static,
    $small.z_vec_static,
)

# This returns kernels
KERNELS
# # Add CUDA kernels
# # NOT SURE HOW THE CUDA syncing works.... putting @CUDA.sync doesn't work before ? 
# # Does the `synchronize()` help things
# if use_cuda
#     KERNELS["large"]["combined_w_mul_cuda"] = @benchmarkable combined_w_mul!(ηϵ, z,
#                                                                              $large.f_w_cuda,
#                                                                              $large.f_cuda,
#                                                                              $large.h_x_cuda,
#                                                                              $large.g_x_cuda,
#                                                                              $large.Q_cuda,
#                                                                              $large.w_0_cuda) setup = (ηϵ = copy(large.ηϵ_cuda); z = copy(large.z_cuda); CUDA.synchronize())
#     KERNELS["medium"]["combined_w_mul_cuda"] = @benchmarkable combined_w_mul!(ηϵ, z,
#                                                                               $medium.f_w_cuda,
#                                                                               $medium.f_cuda,
#                                                                               $medium.h_x_cuda,
#                                                                               $medium.g_x_cuda,
#                                                                               $medium.Q_cuda,
#                                                                               $medium.w_0_cuda) setup = (ηϵ = copy(medium.ηϵ_cuda); z = copy(medium.z_cuda); CUDA.synchronize())
#     KERNELS["small"]["combined_w_mul_cuda"] = @benchmarkable combined_w_mul!(ηϵ, z,
#                                                                              $small.f_w_cuda,
#                                                                              $small.f_cuda,
#                                                                              $small.h_x_cuda,
#                                                                              $small.g_x_cuda,
#                                                                              $small.Q_cuda,
#                                                                              $small.w_0_cuda) setup = (ηϵ = copy(small.ηϵ_cuda); z = copy(small.z_cuda); CUDA.synchronize())
# end

#run(KERNELS)
