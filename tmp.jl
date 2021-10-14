using ChainRulesCore, Zygote

# Will only generate rrule for the `p` argument, and not the p_f!
function f(p, p_f)
    all_p = merge(p, p_f)

    return [all_p.a^2 + all_p.b
            all_p.c^2 + all_p.d^2]
end

function g_1(val)
    p = (a = val[1], b = val[2], d = val[3])
    p_f = (c = 4.0,)
    y = f(p, p_f)
    return sum(y)
end

function g_2(val)
    p = (b = val[1],d = val[2])
    p_f = (a = 5.0, c = 4.0,)
    y = f(p, p_f)
    return sum(y)
end


function ChainRulesCore.rrule(::typeof(f), p, p_f)
    println("using custom rrule")
    val = f(p, p_f)
    # Could go through the keys(p) and fill in details for the derivatives here based on what is in it?
    all_p = merge(p, p_f)
    loadings = [(a = 2 * all_p.a, b = 1, c = 0, d = 0), (a = 0, b = 0, c = 2 * all_p.c, d = 2 * all_p.d)]
    function f_pb(Δsol)
        println("using custom pb")
        Δp = [] # NEED TO FILL IN GIVEN THE APPROPRIATE KEY IN ORDER
        for key in keys(p)
            push!(Δp, loadings[1][key] * Δsol[1] + loadings[2][key] * Δsol[2])
        end
        return NoTangent(), Tangent{typeof(p)}(; zip(keys(p), Δp)...), NoTangent() 
    end
    return val, f_pb
end

# The example down there works
f(nt) = nt.a + 2 * nt.b

function ChainRulesCore.rrule(::typeof(f), nt::NamedTuple)
    y = f(nt)
    # assume some cache passed in here
    tmp = (a = 1, b = 4)
    function f_pullback(dy)
        res_values = []
        for key in keys(nt)
            push!(res_values, tmp[key] * dy)
        end
        # res = NamedTuple{keys(nt)}(res_values)
        # return NoTangent(), Tangent{typeof(nt)}(;res...)
        return NoTangent(), Tangent{typeof(nt)}(; zip(keys(nt), res_values)...)
    end
    # f_pullback(dy) = NoTangent(), Tangent{typeof(nt)}(a=dy, b=3*dy)
    
    return y, f_pullback
end
