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
    solve_quadratic(a, b, c) = (-b + sqrt(b^2 - 4*a*c)) / (2*a)

    LBF_θ((; v0, α_T, p0, Cp), p, T)                                                = T / ( 1 + v0 * α_T * ( p - p0 ) / Cp )
    LBF_temperature_from_pθ((; v0, α_T, p0, Cp), p, θ)                              = θ * ( 1 + v0 * α_T * ( p - p0 ) / Cp )
    LBF_entropy_from_θ((; Cp, T0), θ)                                               = Cp * log( θ  / T0 )
    LBF_θ_from_entropy((; T0, Cp), s)                                               = T0 * exp(s / Cp)
    LBF_specific_volume_from_pθS((; v0, p0, T0, α_T, α_p, α_S), p, θ, S)             = v0 * ( 1 + α_T * ( θ - T0 ) - α_S * (S - S0) - α_p * ( p - p0 ) )
    LBF_sound_speed2((; v0, α_p), T)                                                = (v0 / α_p) * one(T)
    LBF_heat_capacity((; Cp), T)                                                    = Cp * one(T)
    LBF_heat_capacity_s((; Cp), s)                                                  = Cp * one(s)
    LBF_pressure_from_vθS((; p0, T0, v0, α_T, α_p, α_S), v, θ, S)                   = p0 + ( α_T * (θ - T0) - α_S * (S - S0) +  1 - v / v0 ) / α_p  
    LBF_pressure_from_vTS((; Cp, T0, v0, α_T, p0, α_p, α_S), v, T, S)               = p0 + solve_quadratic(v0 * α_T * α_p, α_p * Cp + α_T^2 * T0 * v0 + α_T * ( v - v0 ) + v0 * α_T * α_S * ( S - S0 ), Cp * (α_S * (S - S0) - α_T * (T - T0) + ( v/v0 - 1 )))
    LBF_θ_from_pvS((; Cp, T0, v0, α_T, p0, α_p, α_S), p, v, S)                      = T0 + ( v / v0 - 1 + α_p * (p - p0) + α_S * (S - S0) ) / α_T
    LBF_enthalpy_from_pθS((; Cp, T0, v0, α_T, p0, α_p, α_S, μ0), p, θ, S)           = Cp * (θ - T0) + μ0 * (S - S0) + v0 * (1 + α_T * (θ - T0) - α_S * (S - S0)) * (p - p0) - 0.5 * v0 * α_p * (p - p0)^2
    LBF_gibbs_from_pθS((; Cp, T0, v0, α_T, p0, α_p, α_S, μ0), p, θ, S)              = LBF_enthalpy_from_pθS((; Cp, T0, v0, α_T, p0, α_p), p, θ, S) - LBF_temperature_from_pθ((; v0, α_T, p0, Cp), p, θ) * LBF_entropy_from_θ((; Cp, T0), θ)
    LBF_internal_energy_from_pθS((; Cp, T0, v0, α_T, p0, α_p, α_S, μ0), p, θ, S)    = LBF_enthalpy_from_pθS((; Cp, T0, v0, α_T, p0, α_p, α_S, μ0), p, θ, S) - p * LBF_specific_volume_from_pθS((; v0, p0, T0, α_T, α_p, α_S), p, θ, S)
    LBF_temperature_vsS((; Cp, T0, v0, α_T, p0, α_p), v, s, S)                      = T0 * exp(s / Cp) * ( 1 + (v0 * α_T / (Cp * α_p) ) * ( 1 - v/v0  - α_S * (S - S0) + α_T * T0 * (exp(s / Cp) - 1) ) )
    LBF_exner((; v0, α_T, p0, Cp), p, T)                                            = Cp * T / LBF_θ((; v0, α_T, p0, Cp), p, T)
    LBF_chemical_potential((; μ0, v0, α_S, p0), p)                                  = μ0 - v0 * α_S * ( p - p0 )

    ## all consvar
    # Defined throughout ClimFluids:
    specific_entropy(           fluid::LBF, (p, T, S)::PTSa)    = LBF_entropy_from_θ(fluid, LBF_θ(fluid, p, T))
    specific_enthalpy(          fluid::LBF, (p, T, S)::PTSa)    = LBF_enthalpy_from_pθS(fluid, p, LBF_θ(fluid, p, T), S)
    specific_internal_energy(   fluid::LBF, (p, T, S)::PTSa)    = LBF_internal_energy_from_pθS(fluid, p, LBF_θ(fluid, p, T), S)
    sound_speed2(               fluid::LBF, (p, T, S)::PTSa)    = LBF_sound_speed2(fluid, T)
    potential_temperature(      fluid::LBF, (p, T, S)::PTSa)    = LBF_θ(fluid, p, T)
    potential_enthalpy(         fluid::LBF, (p, T, S)::PTSa)    = LBF_enthalpy_from_pθS(fluid, fluid.p0, LBF_θ(fluid, p, T), S)
    potential_volume(           fluid::LBF, (p, T, S)::PTSa)    = LBF_specific_volume_from_pθS(fluid, fluid.p0, LBF_θ(fluid, p, T), S)
    heat_capacity(              fluid::LBF, (p, T, S)::PTSa)    = LBF_heat_capacity(fluid, T)
    pressure(                   fluid::LBF, (v, T, S)::VTSa)    = LBF_pressure_from_vTS(fluid, v, T, S)
    
    # Defined here only:
    specific_gibbs(             fluid::LBF, (p, T, S)::PTSa)    = LBF_gibbs_from_pθS(fluid, p, LBF_θ(fluid, p, T), S)
    specific_volume(            fluid::LBF, (p, T, S)::PTSa)    = LBF_specific_volume_from_pθS(fluid, p, LBF_θ(fluid, p, T), S)
    specific_entropy(           fluid::LBF, (p, θ, S)::PThSa)   = LBF_entropy_from_θ(fluid, θ)
    chemical_potential(         fluid::LBF, (p, T, S)::PTSa)    = LBF_chemical_potential(fluid, p)
    
    ## consvar = :entropy

    # PTSa
    conservative_variable(      fluid::LBFS, (p, T, S)::PTSa)        = LBF_entropy_from_θ(fluid, LBF_θ(fluid, p, T))
    conjugate_variable(         fluid::LBFS, (p, T, S)::PTSa)        = T
    modified_chemical_potential(fluid::LBFS, (p, T, S)::PTSa)        = LBF_chemical_potential(fluid, p)


    # PConsSa
    temperature(        fluid::LBFS, (p, s, S)::PConsSa)         = LBF_temperature_from_pθ(fluid, p, LBF_θ_from_entropy(fluid, s))
    specific_enthalpy(  fluid::LBFS, (p, s, S)::PConsSa)         = LBF_enthalpy_from_pθS(fluid, p, LBF_θ_from_entropy(fluid, s), S)
    specific_volume(    fluid::LBFS, (p, s, S)::PConsSa)         = LBF_specific_volume_from_pθS(fluid, p, LBF_θ_from_entropy(fluid, s), S)
    heat_capacity(      fluid::LBFS, (p, s, S)::PConsSa)         = LBF_heat_capacity_s(fluid, s)

    function exner_functions(fluid::LBFS, (p, s, S)::PConsSa)    # returns h, v, conjvar = T
        θ = LBF_θ_from_entropy(fluid, s)
        T = LBF_temperature_from_pθ(fluid, p, θ)
        h = LBF_enthalpy_from_pθS(fluid, p, θ, S)
        v = LBF_specific_volume_from_pθS(fluid, p, θ, S)
        return h, v, T
    end
    function volume_functions(fluid::LBFS, (p, s, S)::PConsSa)   # returns v, ∂v/∂p, ∂v/∂s
        θ = LBF_θ_from_entropy(fluid, s)
        v = LBF_specific_volume_from_pθS(fluid, p, θ, S)
        dv_dp = - fluid.α_p * fluid.v0
        dv_ds = θ * fluid.v0 * fluid.α_T / fluid.Cp
        return v, dv_dp, dv_ds
    end
    
    # VConsSa
    temperature(    fluid::LBFS, (v, s)::VConsSa)  = LBF_temperature_vsS(fluid, v, s, S)
    pressure(       fluid::LBFS, (v, s)::VConsSa)  = LBF_pressure_from_vθS(fluid, v, LBF_θ_from_entropy(fluid, s), S)
    
    
    ## consvar = :potential_temperature
    
    # PTSa
    conservative_variable(  fluid::LBFP, (p, T, S)::PTSa)        = LBF_θ(fluid, p, T)
    conjugate_variable(     fluid::LBFP, (p, T, S)::PTSa)        = LBF_exner(fluid, p, T)
    modified_chemical_potential(fluid::LBFS, (p, T, S)::PTSa)    = LBF_chemical_potential(fluid, p)

    # PConsSa
    temperature(        fluid::LBFP, (p, θ, S)::PConsSa)         = LBF_temperature_from_pθ(fluid, p, θ)
    specific_enthalpy(  fluid::LBFP, (p, θ, S)::PConsSa)         = LBF_enthalpy_from_pθS(fluid, p, θ, S)
    specific_volume(    fluid::LBFP, (p, θ, S)::PConsSa)         = LBF_specific_volume_from_pθS(fluid, p, θ, S)
    heat_capacity(      fluid::LBFP, (p, θ, S)::PConsSa)         = LBF_heat_capacity_s(fluid, θ)

    function exner_functions(fluid::LBFP, (p, θ, S)::PConsSa)    # returns h, v, conjvar = Cp0_ct * T / θ
        h = LBF_enthalpy_from_pθS(fluid, p, θ, S)
        v = LBF_specific_volume_from_pθS(fluid, p, θ, S)
        conjvar = fluid.Cp * LBF_temperature_from_pθ(fluid, p, θ) / θ
        return h, v, conjvar
    end
    function volume_functions(fluid::LBFP, (p, θ, S)::PConsSa)   # returns v, ∂v/∂p, ∂v/∂θ
        v = LBF_specific_volume_from_pθS(fluid, p, θ, S)
        dv_dp = - fluid.α_p * fluid.v0
        dv_dθ = fluid.v0 * fluid.α_T
        return v, dv_dp, dv_dθ
    end

    # VCons
    temperature(    fluid::LBFP, (v, θ, S)::VConsSa)    = LBF_temperature_vsS(fluid, v, LBF_entropy_from_θ(fluid, θ), S)
    pressure(       fluid::LBFP, (v, θ, S)::VConsSa)    = LBF_pressure_from_vθS(fluid, v, θ, S)

    ## Fallback: convert to (p, T) if not implemented above
    canonical_state(fluid::LBF, (p, consvar, S)::PConsSa)        = (p, T = temperature(fluid, (; p, consvar, S)), S)
    canonical_state(fluid::LBF, (p, s, S)::PSSa)                 = (p, T = LBF_temperature_from_pθ(fluid, p, LBF_θ_from_entropy(fluid, s)), S)
    canonical_state(fluid::LBF, (p, θ, S)::PThSa)                = (p, T = LBF_temperature_from_pθ(fluid, p, θ), S)
    canonical_state(fluid::LBF, (p, v, S)::PVSa)                 = (p, T = LBF_temperature_from_pθ(fluid, p, LBF_θ_from_pv(fluid, p, v)), S)
    canonical_state(fluid::LBF, (v, T, S)::VTSa)                 = (p = pressure(fluid, (; v, T, S)), T, S)
    canonical_state(fluid::LBF, (v, s, S)::VSSa)                 = canonical_state_LBF_vs(fluid, (v, s, S))
    canonical_state(fluid::LBFS, (v, s, S)::VConsSa)             = canonical_state_LBF_vs(fluid, (v, s, S))
    canonical_state(fluid::LBFP, (v, θ, S)::VConsSa)             = canonical_state_LBF_vθ(fluid, (v, θ, S))

    function canonical_state_LBF_vs(fluid, (v, s, S))
        θ = LBF_θ_from_entropy(fluid, s)
        p = LBF_pressure_from_vθ(fluid, v, θ)
        T = LBF_temperature_from_pθ(fluid, p, θ)
        return (; p, T, S)
    end

    function canonical_state_LBF_vθ(fluid, (v, θ, S))
        p = LBF_pressure_from_vθ(fluid, v, θ)
        T = LBF_temperature_from_pθ(fluid, p, θ)
        return (; p, T, S)
    end

    exner_functions(fluid::LBF, (p, T, S)::PTSa)     = exner_functions(fluid,   (; p, consvar=conservative_variable(fluid, (; p, T, S)), S))
    volume_functions(fluid::LBF, (p, T, S)::PTSa)    = volume_functions(fluid,  (; p, consvar=conservative_variable(fluid, (; p, T, S)), S))

end