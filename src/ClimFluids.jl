module ClimFluids
using MuladdMacro

export AbstractFluid, IdealPerfectGas, CpVarPerfectGas
export temperature,
    specific_entropy, specific_volume, specific_enthalpy, specific_internal_energy
export potential_temperature, conservative_variable, exner_functions
export sound_speed, components

# By default, dot expressions treat user-defined types as iterators, not scalars.
# Derive your types from this abstract type if you want
# them to be treated as scalar in dot expressions.
abstract type ScalarType end
Broadcast.broadcastable(x::ScalarType) = Ref(x)

"""
AbstractFluid{F,N}

Parent type for types describing thermodynamics of a fluid with N components.
"""
abstract type AbstractFluid{F,N} <: ScalarType end
Base.eltype(::AbstractFluid{F}) where {F} = F
@inline components(::AbstractFluid{F,N}) where {F,N} = N

SimpleFluid{F} = AbstractFluid{F,1}
BinaryFluid{F} = AbstractFluid{F,2}

"""
    consvar = conservative_variable(fluid, state)
Returns the conservative variable given a `fluid` and a `state` (named tuple of state variables).
Which conservative variable is chosen is a parameter given to the constructor of `fluid`.
"""
function conservative_variable end

"""
    p = pressure(fluid, state)
Returns pressure given a `fluid` and a `state` (named tuple of state variables).
"""
function pressure end

"""
    T = temperature(fluid, state)
Returns temperature given a `fluid` and a `state` (named tuple of state variables).
"""
function temperature end

"""
    theta = potential_temperature(fluid, state)
Returns potential temperature given a `fluid` and a `state` (named tuple of state variables).
"""
function potential_temperature end

"""
    vpot = potential_volume(fluid, state)
Returns potential volume given a `fluid` and a `state` (named tuple of state variables).
"""
function potential_volume end

"""
    hpot = potential_enthalpy(fluid, state)
Returns potential enthalpy given a `fluid` and a `state` (named tuple of state variables).
"""
function potential_enthalpy end

"""
    v = specific_volume(fluid, state)
Returns specific volume given a `fluid` and a `state` (named tuple of state variables).
"""
function specific_volume end

"""
    s = specific_entropy(fluid, state)
Returns specific entropy given a `fluid` and a `state` (named tuple of state variables).
"""
function specific_entropy end

"""
    h = specific_enthalpy(fluid, state)
Returns specific enthalpy given a `fluid` and a `state` (named tuple of state variables).
"""
function specific_enthalpy end

"""
    h, v, exner, [exner_q] = exner_functions(fluid, state)
Returns specific enthalpy, regarded as a function of pressure, conservative variable and (if applicable) composition,
and its derivatives : specific volume, and Exner-like functions.
"""
function exner_functions end

"""
    e = specific_internal_energy(fluid, state)
Returns specific internal energy given a `fluid` and a `state` (named tuple of state variables).
"""
function specific_internal_energy end

"""
    cs = sound_speed(fluid, (;p,T))
Returns the speed of sound given a `fluid` and a `state` (named tuple of state variables).
Default implementation returns `sqrt(sound_speed2(fluid, state))`.
"""
function sound_speed end

"""
    c2 = sound_speed2(fluid, (;p,T))
Returns the squared speed of sound given a `fluid` and a `state` (named tuple of state variables).
"""
function sound_speed2 end

"""
    cstate = canonical_state(fluid, state)
Converts the named tuple `state` of state variables to the preferred
tuple `cstate`. It is expected that all thermodynamic functions
are implemented for arguments `(fluid, cstate)`. This function is used
internally to provide fallback implementations for thermodynamic functions.
"""
function canonical_state(::F, ::NamedTuple{V}) where {F,V}
    error("""
`canonical_state(::$F, ::NamedTuple{$V})` has been invoked,
indicating that a thermodynamic function was not implemented for fluid $F
and input state variables $V, and a fallback implementation has been called. However
$F does not implement `canonical_state`, on which `ClimFluids`
relies to provide this fallback implementation.
To fix this error, implement :
    ClimFluids.canonical_state(fluid::$F, state::NamedTuple{$V})
""")
end

include("julia/inlineall.jl")
include("julia/state_variables.jl")
include("julia/perfectgas.jl")
include("julia/ideal.jl")
include("julia/lebonnois.jl")
include("julia/binarygas.jl")

@inline pow_fast(x, y) = @fastmath exp(y * log(x))

module ForwardDiff
# implemented only when ForwardDiff is loaded
function exner_functions end
end

#========== for Julia <1.9 ==========#

using PackageExtensionCompat
function __init__()
    @require_extensions
end

end # module
