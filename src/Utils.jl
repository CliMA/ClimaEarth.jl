module Utils

using Printf
using Dates: AbstractTime

maybe_int(t) = isinteger(t) ? Int(t) : t

prettysummary(x, args...) = summary(x)

function prettysummary(f::Function, showmethods=true)
    ft = typeof(f)
    mt = ft.name.mt
    name = mt.name
    n = length(methods(f))
    m = n==1 ? "method" : "methods"
    sname = string(name)
    isself = isdefined(ft.name.module, name) && ft == typeof(getfield(ft.name.module, name))
    ns = (isself || '#' in sname) ? sname : string("(::", ft, ")")
    if showmethods
        return string(ns, " (", "generic function", " with $n $m)")
    else
        return string(ns)
    end
end

prettysummary(x::Int, args...) = string(x)

# This is very important
function prettysummary(nt::NamedTuple, args...)
    n = nfields(nt)

    if n == 0
        return "NamedTuple()"
    else
        str = "("
        for i = 1:n
            f = nt[i]
            str = string(str, fieldname(typeof(nt), i), "=", prettysummary(getfield(nt, i)))
            if n == 1
                str = string(str, ",")
            elseif i < n
                str = string(str, ", ")
            end
        end
    end

    return string(str, ")")
end

function prettykeys(t)
    names = collect(keys(t))
    length(names) == 0 && return "()"
    length(names) == 1 && return string(first(names))
    return string("(", (string(n, ", ") for n in names[1:end-1])..., last(names), ")")
end

# Useful aliases
const second = 1.0
const seconds = second
const minute = 60seconds
const minutes = minute
const hour = 60minutes
const hours = hour
const day = 24hours
const days = day

"""
    prettytime(t, longform=true)

Convert a floating point value `t` representing an amount of time in
SI units of seconds to a human-friendly string with three decimal places.
Depending on the value of `t` the string will be formatted to show `t` in
nanoseconds (ns), microseconds (μs), milliseconds (ms),
seconds, minutes, hours, or days.

With `longform=false`, we use s, m, hrs, and d in place of seconds,
minutes, and hours.
"""
function prettytime(t, longform=true)
    # Modified from: https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/src/trials.jl
    
    # Some shortcuts
    s = longform ? "seconds" : "s" 
    iszero(t) && return "0 $s"
    t < 1e-9 && return @sprintf("%.3e %s", t, s) # yah that's small

    t = maybe_int(t)
    value, units = prettytimeunits(t, longform)

    if isinteger(value)
        return @sprintf("%d %s", value, units)
    else
        return @sprintf("%.3f %s", value, units)
    end
end

function prettytimeunits(t, longform=true)
    t < 1e-9 && return t, "" # just _forget_ picoseconds!
    t < 1e-6 && return t * 1e9, "ns"
    t < 1e-3 && return t * 1e6, "μs"
    t < 1    && return t * 1e3, "ms"
    if t < minute
        value = t
        !longform && return value, "s"
        units = value == 1 ? "second" : "seconds"
        return value, units
    elseif t < hour
        value = maybe_int(t / minute)
        !longform && return value, "m"
        units = value == 1 ? "minute" : "minutes"
        return value, units
    elseif t < day
        value = maybe_int(t / hour)
        units = value == 1 ? (longform ? "hour" : "hr") :
                             (longform ? "hours" : "hrs")
        return value, units
    else
        value = maybe_int(t / day)
        !longform && return value, "d"
        units = value == 1 ? "day" : "days"
        return value, units
    end
end

end # module Utils
