"""
`ClimFluids` provides types and functions to compute
thermodynamic functions of various fluids.

# Using ClimFluids
```
    F = Float32
    N = 1_000_000
    p = F.(1e5 .+ randn(N))
    T = F.(300. .+ randn(N));

    params = (Cp=1000.0, kappa=2/7, p0=1e5, T0=273.15)
    params = (; consvar=:temperature, map(F, params)...)
    fluid = IdealPerfectGas(params)

    s = fluid(:p,:T).specific_entropy.(p,T)
    theta = fluid(:p,:T).conservative_variable.(p,T)
```
In the above example, `fluid` is an object describing an ideal
perfect gas with its physical parameters.
`fluid(:p,:T)` is an object derived from `fluid` which is specialized
for thermodynamics functions taking pressure and temperature as inputs.
`fluid(:p,:T).specific_entropy` is a function-like object computing
the specific entropy of the fluid given `p,T`.
In the named tuple of parameters passed to the constructor `IdealPerfectGas`,
`consvar = :temperature` indicates that potential temperature is chosen as the conservative variable.
As a result the last line computes potential temperature.

See `state_variables(fluid)` for a list of valid inputs to `fluid(...)`
and `state_functions(fluid)` for a list of valid thermodynamic functions.

# Extending ClimFluids
`ClimFluids` is designed to be extensible. To add a new fluid,
define a custom type derived from `SimpleFluid`, `BinaryFluid`, or `AbstractFluid`.
This type should have a constructor accepting a named tuple as single input (see above example).

Thermodynamic functions defined in `ClimFluids` should then be specialized for
this custom type. These functions accept two arguments, one for the fluid and a named
tuple for the input state variables. In the above example, the callable object `fluid(:p,;T).specific_entropy`
calls `ClimFluids.specific_entropy(fluid, (p=p, T=T))`.
Since the same thermodynamic function can be evaluated for different tuples of state
variables, thermodynamic functions must be specialized also for the type of their second argument.
`ClimFluids` provides shortcuts for these types. For example
    `ClimFluids.PT` is the type `NamedTuple{(:p,:T)}`.

Thermodynamic functions may be called in performance-critical inner loops.
It is recommended to make them `@inline` and `@fastmath`.
"""
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

end # module
