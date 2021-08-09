
# Caching utilties
struct ThreadLocalCache{T}
    sl::Threads.SpinLock
    caches::Dict{Int64,T}
end

function (c::ThreadLocalCache)()
    # lock while accessing caches in case deepcopy is slow 
    lock(c.sl) do
        return get!(c.caches, Threads.threadid()) do
            deepcopy(first(values(c.caches))) # copies the first value
        end
    end
end
"""
Thread-local caches, automatically generating new caches when required.

$(SIGNATURES)


# Examples
```julia-repl
julia> caches = ThreadLocalCache(rand(5))

julia> caches()  # get a local cache for the thread 
```
"""
function ThreadLocalCache(cache)
    return ThreadLocalCache(Threads.SpinLock(), Dict(Threads.threadid() => cache))
end

Base.show(io::IO, ::MIME"text/plain", m::ThreadLocalCache) =
    print(io, "ThreadLocalCache with $(m.caches.count) elements\n")

# Cache allocation helper
"""
Helper to allocate threadsafe cache objects consistent with any `AbstractPerturbationModel` for use with `generate_perturbation`.

Initializes with a cache for the current thread.

$(SIGNATURES)


# Examples
```julia-repl
julia> cache = allocate_cache(m);

julia> sol = generate_perturbation(m; cache)  # cache can be used in different threads
```
"""
allocate_cache(m::AbstractPerturbationModel) = ThreadLocalCache(m.allocate_solver_cache(m))

"""
Get a threadsafe cache instance, and determine if the current `p` and `p_f` match it.

NOTE: in the future there could be memoization of multiple caches, but right now there is one per thread.

$(SIGNATURES)


# Examples
```julia-repl
julia> cache = allocate_cache(m);

julia> c, hash_match_all, hash_match_ss, hash_match_perturbation = get_threadsafe_cache(cache, m, p, p_f)
```
"""
function get_threadsafe_cache(
    cache::ThreadLocalCache,
    m::AbstractPerturbationModel,
    p,
    p_f,
    settings = PerturbationSolverSettings(),
)
    # Note later could memoize caches, but unlikely to be helpful here.
    c = cache()  # gets cache associated with threadid(), or allocates as required

    hash_match_all = (hash(p) == c.p_hash) && (hash(p_f) == c.p_f_hash)
    hash_match_ss =
        (hash(get_hash_subset(m.select_p_ss_hash, p)) == c.p_ss_hash) &&
        (hash(get_hash_subset(m.select_p_f_ss_hash, p_f)) == c.p_f_ss_hash)

    hash_match_perturbation =
        (hash(get_hash_subset(m.select_p_perturbation_hash, p)) == c.p_perturbation_hash) &&
        (
            hash(get_hash_subset(m.select_p_f_perturbation_hash, p_f)) ==
            c.p_f_perturbation_hash
        )
    if settings.print_level >= 2 &&
       (hash_match_all || hash_match_ss || hash_match_perturbation)
        println(
            "Cache hit matching: all = $hash_match_all, steady-state = $hash_match_ss, perturbation = $hash_match_perturbation",
        )
        settings.print_level >= 3 && @show (p, p_f)
    end
    return (; c, hash_match_all, hash_match_ss, hash_match_perturbation)
end
