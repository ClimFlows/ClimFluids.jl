"""
    fluid = LinearSimpleFluid(consvar, prec)
    fluid = LinearSimpleFluid(options)
Return an object `fluid` describing the thermodynamics of a single-component simple fluid, linear in potential temperature.
`fluid` can then be used as first argument of thermodynamic functions.
Argument `consvar` defines the conservative variable chosen. Valid values are `:entropy`, `:potential_temperature`.

Arguments can also be passed as a named tuple `options`.
"""
struct LinearSimpleFluid{CV, F} <: SimpleFluid{F}
    p0::F
    T0::F
    Cp::F
    v0::F
    α_T::F
    α_p::F
    function LinearSimpleFluid(consvar::Symbol, (; p0, T0, Cp, v0, α_T, α_p))
        @assert consvar in (:entropy, :potential_temperature)
        new{consvar, typeof(p0)}(p0, T0, Cp, v0, α_T, α_p)
    end
end
LinearSimpleFluid(params) = LinearSimpleFluid(params.consvar, params)

const LSF   = LinearSimpleFluid
const LSFS  = LinearSimpleFluid{:entropy}
const LSFP  = LinearSimpleFluid{:potential_temperature}

@fastmath @muladd @inlineall begin
    ## private functions
    solve_quadratic(a, b, c) = (-b + sqrt(b^2 - 4*a*c)) / (2*a)

    LSF_θ((; v0, α_T, p0, Cp), p, T)                                = T / ( 1 + v0 * α_T * ( p - p0 ) / Cp )
    LSF_temperature_from_pθ((; v0, α_T, p0, Cp), p, θ)              = θ * ( 1 + v0 * α_T * ( p - p0 ) / Cp )
    LSF_entropy_from_θ((; Cp, T0), θ)                               = Cp * log( θ  / T0 )
    LSF_θ_from_entropy((; T0, Cp), s)                               = T0 * exp(s / Cp)
    LSF_specific_volume_from_pθ((; v0, p0, T0, α_T, α_p), p, θ)     = v0 * ( 1 + α_T * ( θ - T0 ) - α_p * ( p - p0 ) )
    LSF_sound_speed2((; v0, α_p), T)                                = (v0 / α_p) * one(T)
    LSF_heat_capacity((; Cp), T)                                    = Cp * one(T)
    LSF_heat_capacity_s((; Cp), s)                                  = Cp * one(s)
    LSF_pressure_from_vθ((; p0, T0, v0, α_T, α_p), v, θ)            = p0 + ( α_T * (θ - T0) +  1 - v / v0 ) / α_p
    LSF_pressure_from_vT((; Cp, T0, v0, α_T, p0, α_p), v, T)        = p0 + solve_quadratic(v0 * α_T * α_p, α_p * Cp + α_T^2 * T0 * v0 + α_T * ( v - v0 ), α_T * Cp * (T0 - T) + Cp * ( v/v0 - 1 ) )
    LSF_θ_from_pv((; Cp, T0, v0, α_T, p0, α_p), p, v)               = T0 + ( v / v0 - 1 + α_p * (p - p0) ) / α_T
    LSF_enthalpy_from_pθ((; Cp, T0, v0, α_T, p0, α_p), p, θ)        = Cp * (θ - T0) + v0 * (1 + α_T * (θ - T0)) * (p - p0) - 0.5 * v0 * α_p * (p - p0)^2
    LSF_gibbs_from_pθ((; Cp, T0, v0, α_T, p0, α_p), p, θ)           = LSF_enthalpy_from_pθ((; Cp, T0, v0, α_T, p0, α_p), p, θ) - LSF_temperature_from_pθ((; v0, α_T, p0, Cp), p, θ) * LSF_entropy_from_θ((; Cp, T0), θ)
    LSF_internal_energy_from_pθ((; Cp, T0, v0, α_T, p0, α_p), p, θ) = LSF_enthalpy_from_pθ((; Cp, T0, v0, α_T, p0, α_p), p, θ) - p * LSF_specific_volume_from_pθ((; v0, p0, T0, α_T, α_p), p, θ)
    LSF_temperature_vs((; Cp, T0, v0, α_T, p0, α_p), v, s)          = T0 *  exp( s / Cp ) * ( 1 + (v0 * α_T / (Cp * α_p) ) * ( 1 - v/v0 + α_T * T0 * (exp(s / Cp) - 1) ) )
    LSF_exner((; v0, α_T, p0, Cp), p, T)                            = Cp * T / LSF_θ((; v0, α_T, p0, Cp), p, T)
    ## all consvar
    # Defined throughout ClimFluids:
    specific_entropy(           fluid::LSF, (p, T)::PT)     = LSF_entropy_from_θ(fluid, LSF_θ(fluid, p, T))
    specific_enthalpy(          fluid::LSF, (p, T)::PT) 	= LSF_enthalpy_from_pθ(fluid, p, LSF_θ(fluid, p, T))
    specific_internal_energy(   fluid::LSF, (p, T)::PT)     = LSF_internal_energy_from_pθ(fluid, p, LSF_θ(fluid, p, T))
    sound_speed2(               fluid::LSF, (p, T)::PT)     = LSF_sound_speed2(fluid, T)
    potential_temperature(      fluid::LSF, (p, T)::PT)     = LSF_θ(fluid, p, T)
    potential_enthalpy(         fluid::LSF, (p, T)::PT)     = LSF_enthalpy_from_pθ(fluid, fluid.p0, LSF_θ(fluid, p, T))
    potential_volume(           fluid::LSF, (p, T)::PT)     = LSF_specific_volume_from_pθ(fluid, fluid.p0, LSF_θ(fluid, p, T))
    heat_capacity(              fluid::LSF, (p, T)::PT)     = LSF_heat_capacity(fluid, T)
    pressure(                   fluid::LSF, (v, T)::VT)     = LSF_pressure_from_vT(fluid, v, T)
    
    # Defined here only:
    specific_gibbs(             fluid::LSF, (p, T)::PT)     = LSF_gibbs_from_pθ(fluid, p, LSF_θ(fluid, p, T))
    specific_volume(            fluid::LSF, (p, T)::PT)     = LSF_specific_volume_from_pθ(fluid, p, LSF_θ(fluid, p, T))
    specific_entropy(           fluid::LSF, (p, θ)::PTh)    = LSF_entropy_from_θ(fluid, θ)

    ## consvar = :entropy

    # PT
    conservative_variable(  fluid::LSFS, (p, T)::PT)        = LSF_entropy_from_θ(fluid, LSF_θ(fluid, p, T))
    conjugate_variable(     fluid::LSFS, (p, T)::PT)        = T

    # PCons
    temperature(        fluid::LSFS, (p, s)::PCons)         = LSF_temperature_from_pθ(fluid, p, LSF_θ_from_entropy(fluid, s))
    specific_enthalpy(  fluid::LSFS, (p, s)::PCons)         = LSF_enthalpy_from_pθ(fluid, p, LSF_θ_from_entropy(fluid, s))
    specific_volume(    fluid::LSFS, (p, s)::PCons)         = LSF_specific_volume_from_pθ(fluid, p, LSF_θ_from_entropy(fluid, s))
    heat_capacity(      fluid::LSFS, (p, s)::PCons)         = LSF_heat_capacity_s(fluid, s)

    function exner_functions(fluid::LSFS, (p, s)::PCons)    # returns h, v, conjvar = T
        θ = LSF_θ_from_entropy(fluid, s)
        T = LSF_temperature_from_pθ(fluid, p, θ)
        h = LSF_enthalpy_from_pθ(fluid, p, θ)
        v = LSF_specific_volume_from_pθ(fluid, p, θ)
        return h, v, T
    end
    function volume_functions(fluid::LSFS, (p, s)::PCons)   # returns v, ∂v/∂p, ∂v/∂s
        θ = LSF_θ_from_entropy(fluid, s)
        v = LSF_specific_volume_from_pθ(fluid, p, θ)
        dv_dp = - fluid.α_p * fluid.v0
        dv_ds = θ * fluid.v0 * fluid.α_T / fluid.Cp
        return v, dv_dp, dv_ds
    end
    
    # VCons
    temperature(    fluid::LSFS, (v, s)::VCons)  = LSF_temperature_vs(fluid, v, s)
    pressure(       fluid::LSFS, (v, s)::VCons)  = LSF_pressure_from_vθ(fluid, v, LSF_θ_from_entropy(fluid, s))
    
    
    ## consvar = :potential_temperature
    
    # PT
    conservative_variable(  fluid::LSFP, (p, T)::PT)        = LSF_θ(fluid, p, T)
    conjugate_variable(     fluid::LSFP, (p, T)::PT)        = LSF_exner(fluid, p, T)

    # PCons
    temperature(        fluid::LSFP, (p, θ)::PCons)         = LSF_temperature_from_pθ(fluid, p, θ)
    specific_enthalpy(  fluid::LSFP, (p, θ)::PCons)         = LSF_enthalpy_from_pθ(fluid, p, θ)
    specific_volume(    fluid::LSFP, (p, θ)::PCons)         = LSF_specific_volume_from_pθ(fluid, p, θ)
    heat_capacity(      fluid::LSFP, (p, θ)::PCons)         = LSF_heat_capacity_s(fluid, θ)

    function exner_functions(fluid::LSFP, (p, θ)::PCons)    # returns h, v, conjvar = Cp0_ct * T / θ
        h = LSF_enthalpy_from_pθ(fluid, p, θ)
        v = LSF_specific_volume_from_pθ(fluid, p, θ)
        conjvar = fluid.Cp * LSF_temperature_from_pθ(fluid, p, θ) / θ
        return h, v, conjvar
    end
    function volume_functions(fluid::LSFP, (p, θ)::PCons)   # returns v, ∂v/∂p, ∂v/∂θ
        v = LSF_specific_volume_from_pθ(fluid, p, θ)
        dv_dp = - fluid.α_p * fluid.v0
        dv_dθ = fluid.v0 * fluid.α_T
        return v, dv_dp, dv_dθ
    end

    # VCons
    temperature(    fluid::LSFP, (v, θ)::VCons)    = LSF_temperature_vs(fluid, v, LSF_entropy_from_θ(fluid, θ))
    pressure(       fluid::LSFP, (v, θ)::VCons)    = LSF_pressure_from_vθ(fluid, v, θ)

    ## Fallback: convert to (p, T) if not implemented above
    canonical_state(fluid::LSF, (p, consvar)::PCons)        = (p, T = temperature(fluid, (; p, consvar)))
    canonical_state(fluid::LSF, (p, s)::PS)                 = (p, T = LSF_temperature_from_pθ(fluid, p, LSF_θ_from_entropy(fluid, s)))
    canonical_state(fluid::LSF, (p, θ)::PTh)                = (p, T = LSF_temperature_from_pθ(fluid, p, θ))
    canonical_state(fluid::LSF, (p, v)::PV)                 = (p, T = LSF_temperature_from_pθ(fluid, p, LSF_θ_from_pv(fluid, p, v)))
    canonical_state(fluid::LSF, (v, T)::VT)                 = (p = pressure(fluid, (; v, T)), T)
    canonical_state(fluid::LSF, (v, s)::VS)                 = canonical_state_LSF_vs(fluid, (v, s))
    canonical_state(fluid::LSFS, (v, s)::VCons)             = canonical_state_LSF_vs(fluid, (v, s))
    canonical_state(fluid::LSFP, (v, θ)::VCons)             = canonical_state_LSF_vθ(fluid, (v, θ))

    function canonical_state_LSF_vs(fluid, (v, s))
        θ = LSF_θ_from_entropy(fluid, s)
        p = LSF_pressure_from_vθ(fluid, v, θ)
        T = LSF_temperature_from_pθ(fluid, p, θ)
        return (; p, T)
    end

    function canonical_state_LSF_vθ(fluid, (v, θ))
        p = LSF_pressure_from_vθ(fluid, v, θ)
        T = LSF_temperature_from_pθ(fluid, p, θ)
        return (; p, T)
    end

    exner_functions(fluid::LSF, (p, T)::PT)     = exner_functions(fluid, (; p, consvar=conservative_variable(fluid, (; p, T))))
    volume_functions(fluid::LSF, (p, T)::PT)    = volume_functions(fluid, (; p, consvar=conservative_variable(fluid, (; p, T))))

end