# IMO this file should be moved downstream to FVGQ20 or 
# somewhere else -- doesn't quite belong in DSSM proper.
# - CSP July 12, 2021

# TODO: remove/modify if something like https://github.com/torfjelde/TuringCallbacks.jl/issues/10 is implemented.
model_name(m::AbstractFirstOrderPerturbationModel{T}) where {T} = "first-$(T)-n$(m.n)"
model_name(m::AbstractSecondOrderPerturbationModel{T}) where {T} = "second-$(T)-n$(m.n)"
function make_turing_callback(
    model;
    stats = nothing,
    comment = "",
    use_tensorboard = true,
    kwargs...,
)  # e.g. stats, include, etc.
    if !use_tensorboard
        return (args...; kwargs...) -> nothing  # i.e. default callback is nothing
    end
    if (comment != "")
        comment = "-" * comment
    end
    if isnothing(stats)
        return TensorBoardCallback(; comment, kwargs...)
    else
        return TensorBoardCallback(stats; comment, kwargs...)
    end
end

function log_turing_results(chain, callback)
    logger = callback.logger  # returns a TBLogger type
    # You can use whatever you want in https://github.com/PhilipVinc/TensorBoardLogger.jl
    p = plot(chain) # or whatever you do for StatsPlots
    with_logger(logger) do
        @info "StatsPlots" p
    end
end
