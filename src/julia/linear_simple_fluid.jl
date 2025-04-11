"""
    fluid = LinearSimpleFluid(consvar, prec)
    fluid = LinearSimpleFluid(options)
Return an object `fluid` describing the thermodynamics of a single-component linear simple fluid.
The parameters are fixed to represent salt-less seawater.
`fluid` can then be used as first argument of thermodynamic functions.
Argument `consvar` defines the conservative variable chosen. Valid values are `:entropy`, 
`:conservative_temperature` corresponding to specific entropy, conservative temperature.

Arguments can also be passed as a named tuple `options`.
"""
struct LinearSimpleFluid{CV, F} <: SimpleFluid{F}
    p0::F
    T0::F
    cp0::F
    v0::F
    cs0::F
    βT::F
    βp::F
    cp0_ct::F

    function LinearSimpleFluid(consvar::Symbol, prec)
        @assert consvar in (:entropy, :conservative_temperature)
        new{consvar, prec}(101325, 273.15, 3986, 9.738e-4, 1490, 1.67e-4, 4.39e-10, 3991.86795711963)
    end
end
LinearSimpleFluid(params) = LinearSimpleFluid(params.consvar, params.prec)

const iter_default = 8     # default number of iterations for Newton method

const LSF   = LinearSimpleFluid
const LSFS  = LinearSimpleFluid{:entropy}
const LSFC  = LinearSimpleFluid{:conservative_temperature}

@fastmath @muladd @inlineall begin
    ## private functions
    LSF_gibbs((; v0, p0, T0, cp0, βp, βT), p, T)            = -cp0 * T * ( log(T * inv(T0)) - 1 ) + v0 * ( p - p0 ) * ( 1 + βT * ( T - T0 ) - 0.5 * βp * ( p - p0 ) )
    LSF_enthalpy((; v0, p0, T0, cp0, βp, βT), p, T)         = cp0 * T + v0 * ( p - p0 ) * ( 1 - βT * T0 - 0.5 * βp * ( p - p0 ) )
    LSF_internal_energy((; v0, p0, T0, cp0, βp, βT), p, T)  = cp0 * T - v0 * ( p0 + βT * ( p * T - p0 * T0 ) - 0.5 * βp * ( p * p - p0 * p0 ) )
    LSF_specific_volume((; v0, p0, T0, βp, βT), p, T)       = v0 * ( 1 + βT * ( T - T0 ) - βp * ( p - p0 ) )    #  gibbs_p
    LSF_entropy((; v0, p0, T0, cp0, βT), p, T)              = cp0 * ( log(T *inv(T0)) ) - v0 * βT *( p - p0 )   # -gibbs_T
    LSF_heat_capacity((; cp0), p, T)                        = cp0 * one(T)
    LSF_gibbs_TT((; cp0), p, T)                             = -cp0 * inv(T)
    LSF_gibbs_pp((; v0, βp), p, T)                          = -v0 * βp * one(T)
    LSF_gibbs_pT((; v0, βT), p, T)                          = v0 *  βT * one(T)
    LSF_sound_speed2((; v0, βT, T0, βp, p0, cp0), p, T)     = v0 * ( 1 + βT * ( T - T0 ) - βp * ( p - p0 ) )^2 * cp0 * inv( βp * cp0 - v0 * βT * βT * T)
    LSF_pt((; v0, p0, cp0, βT), p, T)                       = T * exp( - v0 * βT * ( p - p0 ) * inv( cp0 ) ) 
    LSF_pt_from_ct((; cp0_ct, cp0), Θ)                      = ( cp0_ct * Θ ) * inv( cp0 )
    LSF_entropy_from_pt((; v0, p0, T0, cp0, βT), θ)         = LSF_entropy((; v0, p0, T0, cp0, βT), p0, θ)
    LSF_entropy_from_ct((; cp0_ct, v0, p0, T0, cp0, βT), Θ) = LSF_entropy_from_pt((; v0, p0, T0, cp0, βT), LSF_pt_from_ct((; cp0_ct, cp0), Θ))
    LSF_pressure_vT((; p0, T0, v0, βT, βp), v, T)           = p0 + ( βT * (T - T0) +  1 - v * inv(v0))  * inv(βp)
    LSF_temperature_ps((; cp0, p0, T0, v0, βT), p, s)       = T0 * exp( inv(cp0) * (s + v0 * βT * ( p - p0 ) ) )
    LSF_temperature_pv((; p0, T0, v0, βT, βp), p, v)        = T0 + inv(v0 * βT) * ( v - v0 ) + βp * inv(βT) * ( p - p0 )

    function LSF_temperature_vs(fluid, v, s, T, iters)
        (; cp0, T0, v0, βT, βp) = fluid
        inv_βp = inv(βp)
        for _ in 1:iters        
            f = cp0 * ( log(T * inv(T0)) ) - v0 * βT * βT * inv_βp * (T - T0) - s - βT * inv_βp * ( v0 - v )
            f_T = cp0 * inv(T) - v0 * βT * βT * inv_βp
            err = f * inv(f_T)
            T = T - err
        end
        return T
    end

    function LSF_pressure_vs(fluid, v, s, p, iters)
        (; cp0, p0, T0, v0, βT, βp) = fluid
        inv_βTT0 = inv(βT * T0)
        δp = p - p0
        for _ in 1:iters
            g = cp0 * log( 1 + inv_βTT0 * ( v * inv(v0) - 1 + βp * δp ) ) - v0 * βT * δp - s
            g_T = cp0 * βp * inv( βT * T0 + v * inv(v0) - 1 + βp * δp) - v0 * βT
            err = g * inv(g_T)
            δp = δp - err
        end
        return δp + p0
    end

    ## all consvar
    # Defined throughout ClimFluids:
    specific_entropy(           fluid::LSF, (p, T)::PT)     = LSF_entropy(fluid, p, T)
    specific_enthalpy(          fluid::LSF, (p, T)::PT) 	= LSF_enthalpy(fluid, p, T)
    specific_internal_energy(   fluid::LSF, (p, T)::PT)     = LSF_internal_energy(fluid, p, T) 
    sound_speed2(               fluid::LSF, (p, T)::PT)     = LSF_sound_speed2(fluid, p, T)
    potential_temperature(      fluid::LSF, (p, T)::PT)     = LSF_pt(fluid, p, T)
    potential_enthalpy(         fluid::LSF, (p, T)::PT)     = LSF_enthalpy(fluid, fluid.p0, LSF_pt(fluid, p, T))
    potential_volume(           fluid::LSF, (p, T)::PT)     = LSF_specific_volume(fluid, fluid.p0, LSF_pt(fluid, p, T))
    heat_capacity(              fluid::LSF, (p, T)::PT)     = LSF_heat_capacity(fluid, p, T)
    pressure(                   fluid::LSF, (v, T)::VT)     = LSF_pressure_vT(fluid, v, T)
    
    # Defined here only:
    conservative_temperature(   fluid::LSF, (p, T)::PT)     = LSF_enthalpy(fluid, fluid.p0, LSF_pt(fluid, p, T)) * inv(fluid.cp0_ct)
    specific_gibbs(             fluid::LSF, (p, T)::PT)     = LSF_gibbs(fluid, p, T)
    specific_volume(            fluid::LSF, (p, T)::PT)     = LSF_specific_volume(fluid, p, T)


    ## consvar = :entropy

    # PT
    conservative_variable(  fluid::LSFS, (p, T)::PT)        = LSF_entropy(fluid, p, T)
    conjugate_variable(     fluid::LSFS, (p, T)::PT)        = T

    # PCons
    temperature(        fluid::LSFS, (p, s)::PCons)         = LSF_temperature_ps(fluid, p, s)
    specific_enthalpy(  fluid::LSFS, (p, s)::PCons)         = LSF_enthalpy(fluid, p, LSF_temperature_ps(fluid, p, s))
    specific_volume(    fluid::LSFS, (p, s)::PCons)         = LSF_specific_volume(fluid, p, LSF_temperature_ps(fluid, p, s))
    
    function exner_functions(fluid::LSFS, (p, s)::PCons)    # returns h, v, conjvar = T
        T = LSF_temperature_ps(fluid, p, s)
        h = LSF_enthalpy(fluid, p, T)
        v = LSF_specific_volume(fluid, p, T)
        return h, v, T
    end
    function volume_functions(fluid::LSFS, (p, s)::PCons)   # returns v, ∂v/∂p, ∂v/∂s
        T = LSF_temperature_ps(fluid, p, s)
        v = LSF_specific_volume(fluid, p, T)
        cs2 = LSF_sound_speed2(fluid, p, T)
        dv_dp = - v*v / cs2
        dv_ds = -LSF_gibbs_pT(fluid, p, T) * inv( LSF_gibbs_TT(fluid, p, T) )
        return v, dv_dp, dv_ds
    end
    
    # VCons
    temperature(    fluid::LSFS, (v, s)::VCons; T_ = fluid.T0, iters = iter_default)  = LSF_temperature_vs(fluid, v, s, T_, iters)
    pressure(       fluid::LSFS, (v, s)::VCons; p_ = fluid.p0, iters = iter_default)  = LSF_pressure_vs(fluid, v, s, p_, iters)
    
    
    ## consvar = :conservative_temperature
    
    # PT
    conservative_variable(  fluid::LSFC, (p, T)::PT)        = LSF_enthalpy(fluid, fluid.p0, LSF_pt(fluid, p, T)) * inv(fluid.cp0_ct)
    conjugate_variable(     fluid::LSFC, (p, T)::PT)        = fluid.cp0_ct * T * inv( LSF_pt(fluid, p, T) )

    # PCons
    temperature(        fluid::LSFC, (p, Θ)::PCons)         = LSF_temperature_ps(fluid, p, LSF_entropy_from_ct(fluid, Θ))
    specific_enthalpy(  fluid::LSFC, (p, Θ)::PCons)         = LSF_enthalpy(fluid, p, LSF_temperature_ps(fluid, p, LSF_entropy_from_ct(fluid, Θ)))
    specific_volume(    fluid::LSFC, (p, Θ)::PCons)         = LSF_specific_volume(fluid, p, LSF_temperature_ps(fluid, p, LSF_entropy_from_ct(fluid, Θ)))

    function exner_functions(fluid::LSFC, (p, Θ)::PCons)    # returns h, v, conjvar = cp0_ct * T / θ
        s = LSF_entropy_from_ct(fluid, Θ)
        T = LSF_temperature_ps(fluid, p, s)
        h = LSF_enthalpy(fluid, p, T)
        v = LSF_specific_volume(fluid, p, T)
        conjvar = fluid.cp0_ct * T * inv( LSF_pt(fluid, p, T) )
        return h, v, conjvar
    end
    function volume_functions(fluid::LSFC, (p, Θ)::PCons)   # returns v, ∂v/∂p, ∂v/∂Θ
        s = LSF_entropy_from_ct(fluid, Θ)
        T = LSF_temperature_ps(fluid, p, s)
        v = LSF_specific_volume(fluid, p, T)
        θ = LSF_pt(fluid, p, T)
        cs2 = LSF_sound_speed2(fluid, p, T)
        dv_dp = - v*v / cs2
        # ∂v/∂Θ = ∂v/∂T * ∂T/∂Θ = ∂v/∂T / ( ∂η/∂T * dΘ/dη ) = g_TT / ( -g_TT * θ / cp0_ct)
        dv_dΘ = -LSF_gibbs_pT(fluid, p, T) * fluid.cp0_ct * inv( LSF_gibbs_TT(fluid, p, T) * θ ) 
        return v, dv_dp, dv_dΘ
    end

    # VCons
    temperature(    fluid::LSFC, (v, Θ)::VCons; T_ = fluid.T0, iters = iter_default)    = LSF_temperature_vs(fluid, v, LSF_entropy_from_ct(fluid, Θ), T_, iters)
    pressure(       fluid::LSFC, (v, Θ)::VCons; p_ = fluid.p0, iters = iter_default)    = LSF_pressure_vs(fluid, v, LSF_entropy_from_ct(fluid, Θ), p_, iters)

    ## Fallback: convert to (p, T) if not implemented above
    canonical_state(fluid::LSF, (p, consvar)::PCons)        = (p, T = temperature(fluid, (; p, consvar)))
    canonical_state(fluid::LSF, (p, s)::PS)                 = (p, T = LSF_temperature_ps(fluid, p, s))
    canonical_state(fluid::LSF, (p, v)::PV)                 = (p, T = LSF_temperature_pv(fluid, p, v))
    canonical_state(fluid::LSF, (v, T)::VT)                 = (p = pressure(fluid, (; v, T)), T)
    canonical_state(fluid::LSF, (v, s)::VS)                 = canonical_state_LSF_vs(fluid, (v, s))
    canonical_state(fluid::LSFS, (v, s)::VCons)             = canonical_state_LSF_vs(fluid, (v, s))
    canonical_state(fluid::LSFC, (v, Θ)::VCons)             = canonical_state_LSF_vΘ(fluid, (v, Θ))
    
    function canonical_state_LSF_vs(fluid, (v, s))
        T = LSF_temperature_vs(fluid, v, s, fluid.T0, iter_default)
        p = LSF_pressure_vs(fluid, v, s, fluid.p0, iter_default)
        return (; p, T)
    end

    function canonical_state_LSF_vΘ(fluid, (v, Θ))
        s = LSF_entropy_from_ct(fluid, Θ)
        T = LSF_temperature_vs(fluid, v, s, fluid.T0, iter_default)
        p = LSF_pressure_vs(fluid, v, s, fluid.p0, iter_default)
        return (; p, T)
    end

    exner_functions(fluid::LSF, (p, T)::PT)     = exner_functions(fluid, (; p, consvar=conservative_variable(fluid, (; p, T))))
    volume_functions(fluid::LSF, (p, T)::PT)    = volume_functions(fluid, (; p, consvar=conservative_variable(fluid, (; p, T))))

end