"""
    fluid = LinearBinaryFluid(consvar, prec)
    fluid = LinearBinaryFluid(options)
Return an object `fluid` describing the thermodynamics of a binary fluid, linear in potential temperature and salinity.
`fluid` can then be used as first argument of thermodynamic functions.
Argument `consvar` defines the conservative variable chosen. Valid values are `:entropy`, `:potential_temperature`.

Arguments can also be passed as a named tuple `options`.
"""
struct LinearBinaryFluid{CV, F} <: BinaryFluid{F}
    p0::F
    T0::F
    S0::F
    Cp::F
    μ0::F
    v0::F
    α_T::F
    α_S::F
    α_p::F
    function LinearBinaryFluid(consvar::Symbol, (; p0, T0, S0, Cp, μ0, v0, α_T, α_S, α_p))
        @assert consvar in (:entropy, :potential_temperature)
        new{consvar, typeof(p0)}(p0, T0, S0, Cp, μ0, v0, α_T, α_S, α_p)
    end
end
LinearBinaryFluid(params) = LinearBinaryFluid(params.consvar, params)

const LBF   = LinearBinaryFluid
const LBFS  = LinearBinaryFluid{:entropy}
const LBFP  = LinearBinaryFluid{:potential_temperature}

@fastmath @muladd @inlineall begin
    ## private functions
    # solve_quadratic(a, b, c) = (-b + sqrt(b^2 - 4*a*c)) / (2*a)

    LBF_θ((; v0, α_T, p0, Cp), p, T)                                                = T / ( 1 + v0 * α_T * ( p - p0 ) / Cp )
    LBF_temperature_from_pθ((; v0, α_T, p0, Cp), p, θ)                              = θ * ( 1 + v0 * α_T * ( p - p0 ) / Cp )
    LBF_entropy_from_θ((; Cp, T0), θ)                                               = Cp * log( θ  / T0 )
    LBF_θ_from_entropy((; T0, Cp), s)                                               = T0 * exp(s / Cp)
    LBF_specific_volume_from_pθS((; v0, p0, T0, α_T, α_p, α_S), p, θ, q)             = v0 * ( 1 + α_T * ( θ - T0 ) - α_S * (q - S0) - α_p * ( p - p0 ) )
    LBF_sound_speed2((; v0, α_p), T)                                                = (v0 / α_p) * one(T)
    LBF_heat_capacity((; Cp), T)                                                    = Cp * one(T)
    LBF_heat_capacity_s((; Cp), s)                                                  = Cp * one(s)
    LBF_pressure_from_vθS((; p0, T0, v0, α_T, α_p, α_S), v, θ, q)                   = p0 + ( α_T * (θ - T0) - α_S * (q - S0) +  1 - v / v0 ) / α_p  
    LBF_pressure_from_vTS((; Cp, T0, v0, α_T, p0, α_p, α_S), v, T, q)               = p0 + solve_quadratic(v0 * α_T * α_p, α_p * Cp + α_T^2 * T0 * v0 + α_T * ( v - v0 ) + v0 * α_T * α_S * ( q - S0 ), Cp * (α_S * (q - S0) - α_T * (T - T0) + ( v/v0 - 1 )))
    LBF_θ_from_pvS((; Cp, T0, v0, α_T, p0, α_p, α_S), p, v, q)                      = T0 + ( v / v0 - 1 + α_p * (p - p0) + α_S * (q - S0) ) / α_T
    LBF_enthalpy_from_pθS((; Cp, T0, S0, v0, α_T, p0, α_p, α_S, μ0), p, θ, q)           = Cp * (θ - T0) + μ0 * (q - S0) + v0 * (1 + α_T * (θ - T0) - α_S * (q - S0)) * (p - p0) - 0.5 * v0 * α_p * (p - p0)^2
    LBF_gibbs_from_pθS((; Cp, T0, v0, α_T, p0, α_p, α_S, μ0), p, θ, q)              = LBF_enthalpy_from_pθS((; Cp, T0, v0, α_T, p0, α_p), p, θ, q) - LBF_temperature_from_pθ((; v0, α_T, p0, Cp), p, θ) * LBF_entropy_from_θ((; Cp, T0), θ)
    LBF_internal_energy_from_pθS((; Cp, T0, v0, α_T, p0, α_p, α_S, μ0), p, θ, q)    = LBF_enthalpy_from_pθS((; Cp, T0, v0, α_T, p0, α_p, α_S, μ0), p, θ, q) - p * LBF_specific_volume_from_pθS((; v0, p0, T0, α_T, α_p, α_S), p, θ, q)
    LBF_temperature_vsS((; Cp, T0, v0, α_T, p0, α_p), v, s, q)                      = T0 * exp(s / Cp) * ( 1 + (v0 * α_T / (Cp * α_p) ) * ( 1 - v/v0  - α_S * (q - S0) + α_T * T0 * (exp(s / Cp) - 1) ) )
    LBF_exner((; v0, α_T, p0, Cp), p, T)                                            = Cp * T / LBF_θ((; v0, α_T, p0, Cp), p, T)
    LBF_chemical_potential((; μ0, v0, α_S, p0), p, q)                               = μ0 * one(q) - v0 * α_S * ( p - p0 * one(q) )
    LBF_ds_dq(s) = zero(s)
    LBF_dθ_ds((; Cp), θ) = θ * inv(Cp)
    LBF_dθ_dq(θ) = zero(θ)

    ## all consvar
    # Defined throughout ClimFluids:
    specific_entropy(           fluid::LBF, (p, T, q)::PTQ)    = LBF_entropy_from_θ(fluid, LBF_θ(fluid, p, T))
    specific_enthalpy(          fluid::LBF, (p, T, q)::PTQ)    = LBF_enthalpy_from_pθS(fluid, p, LBF_θ(fluid, p, T), q)
    specific_internal_energy(   fluid::LBF, (p, T, q)::PTQ)    = LBF_internal_energy_from_pθS(fluid, p, LBF_θ(fluid, p, T), q)
    sound_speed2(               fluid::LBF, (p, T, q)::PTQ)    = LBF_sound_speed2(fluid, T)
    potential_temperature(      fluid::LBF, (p, T, q)::PTQ)    = LBF_θ(fluid, p, T)
    potential_enthalpy(         fluid::LBF, (p, T, q)::PTQ)    = LBF_enthalpy_from_pθS(fluid, fluid.p0, LBF_θ(fluid, p, T), q)
    potential_volume(           fluid::LBF, (p, T, q)::PTQ)    = LBF_specific_volume_from_pθS(fluid, fluid.p0, LBF_θ(fluid, p, T), q)
    heat_capacity(              fluid::LBF, (p, T, q)::PTQ)    = LBF_heat_capacity(fluid, T)
    pressure(                   fluid::LBF, (v, T, q)::VTQ)    = LBF_pressure_from_vTS(fluid, v, T, q)
    
    # Defined here only:
    specific_gibbs(             fluid::LBF, (p, T, q)::PTQ)    = LBF_gibbs_from_pθS(fluid, p, LBF_θ(fluid, p, T), q)
    specific_volume(            fluid::LBF, (p, T, q)::PTQ)    = LBF_specific_volume_from_pθS(fluid, p, LBF_θ(fluid, p, T), q)
    specific_entropy(           fluid::LBF, (p, θ, q)::PThQ)   = LBF_entropy_from_θ(fluid, θ)
    chemical_potential(         fluid::LBF, (p, T, q)::PTQ)    = LBF_chemical_potential(fluid, p, q)
    ds_dq(                      fluid::LBF, (p, s, q)::PSQ)    = LBF_ds_dq(s)
    ## consvar = :entropy

    # PTQ
    conservative_variable(      fluid::LBFS, (p, T, q)::PTQ)        = LBF_entropy_from_θ(fluid, LBF_θ(fluid, p, T))
    conjugate_variable(         fluid::LBFS, (p, T, q)::PTQ)        = T
    modified_chemical_potential(fluid::LBFS, (p, T, q)::PTQ)        = LBF_chemical_potential(fluid, p, q)

    # PConsQ
    temperature(        fluid::LBFS, (p, s, q)::PConsQ)             = LBF_temperature_from_pθ(fluid, p, LBF_θ_from_entropy(fluid, s))
    specific_enthalpy(  fluid::LBFS, (p, s, q)::PConsQ)             = LBF_enthalpy_from_pθS(fluid, p, LBF_θ_from_entropy(fluid, s), q)
    specific_volume(    fluid::LBFS, (p, s, q)::PConsQ)             = LBF_specific_volume_from_pθS(fluid, p, LBF_θ_from_entropy(fluid, s), q)
    heat_capacity(      fluid::LBFS, (p, s, q)::PConsQ)             = LBF_heat_capacity_s(fluid, s)
    modified_chemical_potential(fluid::LBFS, (p, s, q)::PConsQ)     = LBF_chemical_potential(fluid, p, q)
    dcons_dq(                   fluid::LBFS, (p, s, q)::PConsQ)     = LBF_ds_dq(s)

    function exner_functions(fluid::LBFS, (p, s, q)::PConsQ)    # returns h, v, conjvar = T
        θ = LBF_θ_from_entropy(fluid, s)
        T = LBF_temperature_from_pθ(fluid, p, θ)
        h = LBF_enthalpy_from_pθS(fluid, p, θ, q)
        v = LBF_specific_volume_from_pθS(fluid, p, θ, q)
        return h, v, T
    end
    function volume_functions(fluid::LBFS, (p, s, q)::PConsQ)   # returns v, ∂v/∂p, ∂v/∂s
        θ = LBF_θ_from_entropy(fluid, s)
        v = LBF_specific_volume_from_pθS(fluid, p, θ, q)
        dv_dp = - fluid.v0 * fluid.α_p
        dv_ds = θ * fluid.v0 * fluid.α_T / fluid.Cp
        dv_dS = - fluid.v0 * fluid.α_S
        return v, dv_dp, dv_ds, dv_dS
    end
    
    # VConsQ
    temperature(    fluid::LBFS, (v, s)::VConsQ)  = LBF_temperature_vsS(fluid, v, s, q)
    pressure(       fluid::LBFS, (v, s)::VConsQ)  = LBF_pressure_from_vθS(fluid, v, LBF_θ_from_entropy(fluid, s), q)
    
    
    ## consvar = :potential_temperature
    
    # PTQ
    conservative_variable(  fluid::LBFP, (p, T, q)::PTQ)       = LBF_θ(fluid, p, T)
    conjugate_variable(     fluid::LBFP, (p, T, q)::PTQ)       = LBF_exner(fluid, p, T)
    modified_chemical_potential(fluid::LBFP, (p, T, q)::PTQ)   = LBF_chemical_potential(fluid, p, q)

    # PConsQ
    temperature(        fluid::LBFP, (p, θ, q)::PConsQ)        = LBF_temperature_from_pθ(fluid, p, θ)
    specific_enthalpy(  fluid::LBFP, (p, θ, q)::PConsQ)        = LBF_enthalpy_from_pθS(fluid, p, θ, q)
    specific_volume(    fluid::LBFP, (p, θ, q)::PConsQ)        = LBF_specific_volume_from_pθS(fluid, p, θ, q)
    heat_capacity(      fluid::LBFP, (p, θ, q)::PConsQ)        = LBF_heat_capacity_s(fluid, θ)
    modified_chemical_potential(fluid::LBFP, (p, θ, q)::PConsQ)   = LBF_chemical_potential(fluid, p, q)
    dcons_dq(               fluid::LBFP, (p, θ, q)::PConsQ)   = LBF_dθ_ds(fluid, θ) * LBF_ds_dq(LBF_entropy_from_θ(fluid, θ)) + LBF_dθ_dq(θ)

    function exner_functions(fluid::LBFP, (p, θ, q)::PConsQ)    # returns h, v, conjvar = Cp0_ct * T / θ
        h = LBF_enthalpy_from_pθS(fluid, p, θ, q)
        v = LBF_specific_volume_from_pθS(fluid, p, θ, q)
        conjvar = fluid.Cp * LBF_temperature_from_pθ(fluid, p, θ) / θ
        return h, v, conjvar
    end
    function volume_functions(fluid::LBFP, (p, θ, q)::PConsQ)   # returns v, ∂v/∂p, ∂v/∂θ
        v = LBF_specific_volume_from_pθS(fluid, p, θ, q)
        dv_dp = -fluid.v0 * fluid.α_p
        dv_dθ = fluid.v0 * fluid.α_T
        dv_dS = fluid.v0 * fluid.α_S
        return v, dv_dp, dv_dθ, dv_dS
    end

    # VCons
    temperature(    fluid::LBFP, (v, θ, q)::VConsQ)    = LBF_temperature_vsS(fluid, v, LBF_entropy_from_θ(fluid, θ), q)
    pressure(       fluid::LBFP, (v, θ, q)::VConsQ)    = LBF_pressure_from_vθS(fluid, v, θ, q)

    ## Fallback: convert to (p, T) if not implemented above
    canonical_state(fluid::LBF, (p, consvar, q)::PConsQ)        = (p, T = temperature(fluid, (; p, consvar, q)), q)
    canonical_state(fluid::LBF, (p, s, q)::PSQ)                 = (p, T = LBF_temperature_from_pθ(fluid, p, LBF_θ_from_entropy(fluid, s)), q)
    canonical_state(fluid::LBF, (p, θ, q)::PThQ)                = (p, T = LBF_temperature_from_pθ(fluid, p, θ), q)
    canonical_state(fluid::LBF, (p, v, q)::PVQ)                 = (p, T = LBF_temperature_from_pθ(fluid, p, LBF_θ_from_pv(fluid, p, v)), q)
    canonical_state(fluid::LBF, (v, T, q)::VTQ)                 = (p = pressure(fluid, (; v, T, q)), T, q)
    canonical_state(fluid::LBF, (v, s, q)::VSQ)                 = canonical_state_LBF_vs(fluid, (v, s, q))
    canonical_state(fluid::LBFS, (v, s, q)::VConsQ)             = canonical_state_LBF_vs(fluid, (v, s, q))
    canonical_state(fluid::LBFP, (v, θ, q)::VConsQ)             = canonical_state_LBF_vθ(fluid, (v, θ, q))

    function canonical_state_LBF_vs(fluid, (v, s, q))
        θ = LBF_θ_from_entropy(fluid, s)
        p = LBF_pressure_from_vθ(fluid, v, θ)
        T = LBF_temperature_from_pθ(fluid, p, θ)
        return (; p, T, q)
    end

    function canonical_state_LBF_vθ(fluid, (v, θ, q))
        p = LBF_pressure_from_vθ(fluid, v, θ)
        T = LBF_temperature_from_pθ(fluid, p, θ)
        return (; p, T, q)
    end

    exner_functions(fluid::LBF, (p, T, q)::PTQ)     = exner_functions(fluid,   (; p, consvar=conservative_variable(fluid, (; p, T, q)), q))
    volume_functions(fluid::LBF, (p, T, q)::PTQ)    = volume_functions(fluid,  (; p, consvar=conservative_variable(fluid, (; p, T, q)), q))

end