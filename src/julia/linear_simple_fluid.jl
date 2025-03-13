"""
Linear Simple Fluid EOS
"""
struct LinearSimpleFluid{CV, F} <: SimpleFluid{F}
    p0::F
    T0::F
    cp0::F
    ρ0::F
    α0::F
    cs0::F
    βT::F
    βp::F
    cp0_ct::F
    function LinearSimpleFluid(consvar::Symbol, prec)
        new{consvar, prec}(101325, 273.15, 3986, 1027.0, 9.738e-4, 1490, 1.67e-4, 4.39e-10, 3991.86795711963)
    end
end
LinearSimpleFluid(params) = LinearSimpleFluid(params.consvar, params.prec)

LSF = LinearSimpleFluid
LSFE = LinearSimpleFluid{:entropy}
LSFC = LinearSimpleFluid{:conservative_temperature}

VCons = NamedTuple{(:v, :consvar)}
PCons = NamedTuple{(:p, :consvar)}
PT = NamedTuple{(:p, :T)}
VT = NamedTuple{(:v, :T)}

@fastmath @muladd @inlineall begin
    ## private functions
    LSF_gibbs((; α0, p0, T0, cp0, βp, βT), p, T)            = -cp0 * T * ( log(T * inv(T0)) - 1 ) + α0 * ( p - p0 ) * ( 1 + βT * ( T - T0 ) - 0.5 * βp * ( p - p0 ) )
    LSF_specific_volume((; α0, p0, T0, βp, βT), p, T)       = α0 * ( 1 + βT * ( T - T0 ) - βp * ( p - p0 ) )
    LSF_specific_entropy((; α0, p0, T0, cp0, βT), p, T)     = cp0 * ( log(T *inv(T0)) ) - α0 * βT *( p - p0 ) 
    LSF_specific_heat_capacity((; cp0), p, T)               = cp0 * one(T)
    LSF_gibbs_TT((; cp0), p, T)                             = -cp0 * inv(T)
    LSF_gibbs_pp((; α0, βp), p, T)                          = - α0 * βp * one(T)
    LSF_gibbs_pT((; α0, βT), p, T)                          = α0 *  βT * one(T)
    LSF_specific_enthalpy((; α0, p0, T0, cp0, βp, βT), p, T) = cp0 * T + α0 * ( p - p0 ) * ( 1 - βT * T0 - 0.5 * βp * ( p - p0 ) )
    LSF_potential_temperature((; α0, p0, cp0, βT), p, T)    = T * exp( - α0 * βT * ( p - p0 ) * inv( cp0 ) ) 
    LSF_pt_from_ct((; cp0_ct, cp0), Θ)                      = ( cp0_ct * Θ ) * inv( cp0 )
    LSF_entropy_from_pt((; α0, p0, T0, cp0, βT), θ)         = LSF_specific_entropy((; α0, p0, T0, cp0, βT), p0, θ)
    LSF_entropy_from_ct((; cp0_ct, g0, η0, α0, p0, T0, cp0, βT, η0), Θ) = LSF_entropy_from_pt((; α0, p0, T0, cp0, βT), LSF_pt_from_ct((; cp0_ct, cp0), Θ))
    
    LSF_sound_speed2((; α0, βT, T0, βp, p0, cp0), p, T) = α0 * ( 1 + βT * ( T - T0 ) - βp * ( p - p0 ) )^2 * cp0 * inv( βp * cp0 - α0 * βT * βT * T)
    
    # pressure and temperature from (v, η)
    function LSF_pressure_temperature(fluid, v, η, T, iters)
        (; cp0, p0, T0, α0, βT, βp) = fluid
        inv_βp = inv(βp)
        for i in 1:iters        
            f = cp0 * ( log(T * inv(T0)) )- α0 * βT * βT * inv_βp * (T - T0) - η - βT * inv_βp * ( α0 - v )
            f_T = cp0 * inv(T) - α0 * βT * βT * inv_βp
            err = f * inv(f_T)
            T = T - err
            if abs(err) < 1e-12
                break
            elseif i == iters
                @info err
                error("Newton iteration did not converge")
            end
        end
        p = p0 + inv_βp * ( 1 + βT * ( T - T0 ) - v * inv(α0) )
        return p, T
    end

    function LSF_temperature_v_η(fluid, v, η, T, iters)
        (; cp0, T0, α0, βT, βp) = fluid
        inv_βp = inv(βp)
        for _ in 1:iters        
            f = cp0 * ( log(T * inv(T0)) )- α0 * βT * βT * inv_βp * (T - T0) - η - βT * inv_βp * ( α0 - v )
            f_T = cp0 * inv(T) - α0 * βT * βT * inv_βp
            err = f * inv(f_T)
            T = T - err
        end
        return T
    end
    # function temp_lambert(fluid, v, η)
    #     (; α0, βT, βp, cp0, T0) = fluid
    #     pref = - α0 * βT * βT * inv(βp) * inv(cp0)
    #     T0_ = T0 * pref
    #     inv_βp, inv_cp0 = inv(βp), inv(cp0)
    #     T_ = lambertw( T0_ * exp(T0_ + βT * inv_βp * inv_cp0 * (α0 - v) + η * inv_cp0))
    #     return T_ * inv(pref)
    # end

    # temperature from (p, η)
    LSF_temperature((; cp0, p0, T0, α0, βT), p, η)  = T0 * exp( inv(cp0) * (η + α0 * βT * ( p - p0 ) ) )

    # pressure from (v, T)
    LSF_pressure((; p0, T0, α0, βT, βp), v, T)      = p0 + inv( βp ) * ( 1 + βT * ( T - T0 ) ) - v * inv( α0 * βp ) 

    ## all consvar
    specific_enthalpy(fluid::LSF, (p, T)::PT)   = LSF_specific_enthalpy(fluid, p, T)
    gibbs_function(fluid::LSF, (p, T)::PT)      = LSF_gibbs(fluid, p, T)
    specific_entropy(fluid::LSF, (p, T)::PT)    = LSF_specific_entropy(fluid, p, T)
    specific_volume(fluid::LSF, (p, T)::PT)     = LSF_specific_volume(fluid, p, T)
    sound_speed(fluid::LSF, (p, T)::PT)         = sqrt(LSF_sound_speed2(fluid, p, T))
    sound_speed2(fluid::LSF, (p, T)::PT)        = LSF_sound_speed2(fluid, p, T)
    isobaric_heat_capacity(fluid::LSF, (p, T)::PT)   = LSF_specific_heat_capacity(fluid, p, T)
    pressure(fluid::LSF, (v, T)::VT)            = LSF_pressure(fluid, v, T)

    ## consvar = :entropy

    # PCons
    temperature(fluid::LSFE, (p, η)::PCons)         = LSF_temperature(fluid, p, η)
    specific_enthalpy(fluid::LSFE, (p, η)::PCons)   = LSF_specific_enthalpy(fluid, p, LSF_temperature(fluid, p, η))
    specific_volume(fluid::LSFE, (p, η)::PCons)     = LSF_specific_volume(fluid, p, LSF_temperature(fluid, p, η))
    dT_consvar(fluid::LSFE, (p, η)::PCons)          = LSF_dT_η(fluid, p, η)

    # VCons
    temperature(fluid::LSFE, (v, η)::VCons; T_ = fluid.T0, iters = 10)          = LSF_temperature_v_η(fluid, v, η, T_, iters)
    # pressure_temperature(fluid::LSFE, (v, η)::VCons; T_ = fluid.T0, iters = 3)  = LSF_pressure_temperature(fluid, v, η, T_, iters)
    pressure(fluid::LSFE, (v, η)::VCons; T_ = fluid.T0, iters = 10)              = LSF_pressure(fluid, v, LSF_temperature_v_η(fluid, v, η, T_, iters))
    
    # PT
    conjugate_variable(fluid::LSFE, (p, T)::PT)     = T
    conservative_variable(fluid::LSFE, (p, T)::PT)  = LSF_specific_entropy(fluid, p, T)
    gibbs_pT(fluid::LSFE, (p, T)::PT)               = LSF_gibbs_pT(fluid, p, T)
    gibbs_TT(fluid::LSFE, (p, T)::PT)               = LSF_gibbs_TT(fluid, p, T)
    
    ## consvar = :conservative_temperature

    # PCons
    temperature(fluid::LSFC, (p, Θ)::PCons)         = LSF_temperature(fluid, p, LSF_entropy_from_ct(fluid, Θ))
    specific_enthalpy(fluid::LSFC, (p, Θ)::PCons)   = LSF_specific_enthalpy(fluid, p, LSF_temperature(fluid, p, LSF_entropy_from_ct(fluid, Θ)))
    specific_volume(fluid::LSFC, (p, Θ)::PCons)     = LSF_specific_volume(fluid, p, LSF_temperature(fluid, p, LSF_entropy_from_ct(fluid, Θ)))

    # VCons
    temperature(fluid::LSFC, (v, Θ)::VCons; T_ = fluid.T0, iters = 10)   = LSF_temperature_v_η(fluid, v, LSF_entropy_from_ct(fluid, Θ), T_, iters)
    # pressure_temperature(fluid::LSFC, (v, Θ)::VCons; T_ = fluid.T0, iters = 10) = LSF_pressure_temperature(fluid, v, LSF_entropy_from_ct(fluid, Θ), T_, iters)
    pressure(fluid::LSFC, (v, Θ)::VCons; T_ = fluid.T0, iters = 10) = LSF_pressure(fluid, v, LSF_temperature_v_η(fluid, v, LSF_entropy_from_ct(fluid, Θ), T_, iters))

    # PT
    conjugate_variable(fluid::LSFC, (p, T)::PT)     = fluid.cp0_ct * T * inv( LSF_potential_temperature(fluid, p, T) )
    conservative_variable(fluid::LSFC, (p, T)::PT)  = LSF_specific_enthalpy(fluid, fluid.p0, LSF_potential_temperature(fluid, p, T)) * inv(fluid.cp0_ct)
    
end