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
    val = f(p, p_f)
    # Could go through the keys(p) and fill in details for the derivatives here based on what is in it?

    function f_pb(Δsol)
        Δp = # NEED TO FILL IN GIVEN THE APPROPRIATE KEY IN ORDER
        return NoTangent(), Δp, NoTangent() 
    end
end

# The example down there works
f(nt) = nt.a + 2 * nt.b

function ChainRulesCore.rrule(::typeof(f), nt::NamedTuple)
    y = f(nt)
    f_pullback(dy) = NoTangent(), Tangent{typeof(nt)}(a=dy, b=3*dy)
    
    return y, f_pullback
end