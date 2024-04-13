# ClimFluids

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ClimFlows.github.io/ClimFluids.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ClimFlows.github.io/ClimFluids.jl/dev/)
[![Build Status](https://github.com/ClimFlows/ClimFluids.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ClimFlows/ClimFluids.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ClimFlows/ClimFluids.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ClimFlows/ClimFluids.jl)

`ClimFluids` provides types and functions to compute
thermodynamic functions of various fluids of interest for climate modelling.

## Installing
`ClimFluids` is registered with the ClimFluids registry. See [instructions there](https://github.com/ClimFlows/JuliaRegistry), then:
```julia
]add ClimFluids
```

## Using ClimFluids
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

## Extending ClimFluids
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
