
# Caching utilties

# Cache allocation helper


allocate_cache(m::AbstractPerturbationModel) = m.allocate_solver_cache(m)

# Functions to map to hash subsets
get_hash_subset(A::AbstractMatrix, p) = A * p
get_hash_subset(::UniformScaling, p) = p


# Removing cache checking support for now.
# Would be helpful with gibbs blocks down the road

# function check_cache(
#     cache,
#     m,
#     p,
#     p_f,
#     settings = PerturbationSolverSettings(),
# )
#     hash_match_all = (hash(p) == cache.p_hash) && (hash(p_f) == cache.p_f_hash)
#     hash_match_ss =
#         (hash(get_hash_subset(m.select_p_ss_hash, p)) == cache.p_ss_hash) &&
#         (hash(get_hash_subset(m.select_p_f_ss_hash, p_f)) == cache.p_f_ss_hash)

#     hash_match_perturbation =
#         (hash(get_hash_subset(m.select_p_perturbation_hash, p)) == cache.p_perturbation_hash) &&
#         (
#             hash(get_hash_subset(m.select_p_f_perturbation_hash, p_f)) ==
#             cache.p_f_perturbation_hash
#         )
#     if settings.use_solution_cache && settings.print_level >= 2 &&
#        (hash_match_all || hash_match_ss || hash_match_perturbation)       
#         println(
#             "Cache hit matching: all = $hash_match_all, steady-state = $hash_match_ss, perturbation = $hash_match_perturbation",
#         )
#         settings.print_level >= 3 && @show (p, p_f)
#     end
#     return (; hash_match_all, hash_match_ss, hash_match_perturbation)
# end

# function reset_cache_hashes(c)
#     c.p_hash = zero(UInt64)
#     c.p_f_hash = zero(UInt64)
#     c.p_ss_hash = zero(UInt64)
#     c.p_f_ss_hash = zero(UInt64)
#     c.p_perturbation_hash = zero(UInt64)
#     c.p_f_perturbation_hash = zero(UInt64)
# end
