"""
    fluid = NonlinearBinaryFluid(consvar, prec)
    fluid = NonlinearBinaryFluid(options)
Return an object `fluid` describing the thermodynamics of a nonlinear binary fluid.
The fluid has a constant heat capacity, but includes cabelling and thermobaric terms.
`fluid` can then be used as first argument of thermodynamic functions.
Argument `consvar` defines the conservative variable chosen. Valid values are `:entropy`, `:potential_temperature`.

Arguments can also be passed as a named tuple `options`.
"""
struct NonlinearBinaryFluid{CV, F} <: BinaryFluid{F}
    p0::F
    T0::F
    S0::F
    Cp::F
    μ0::F
    v0::F
    α_T::F
    α_S::F
    α_p::F
    α_TT::F
    γ::F
    Cp0::F
    g0::F
    η0::F
    function NonlinearBinaryFluid(consvar::Symbol, (; p0, T0, S0, Cp, μ0, v0, α_T, α_S, α_p, α_TT, γ))
        @assert consvar in (:entropy, :potential_temperature)
        new{consvar, typeof(p0)}(p0, T0, S0, Cp, μ0, v0, α_T, α_S, α_p, α_TT, γ, 3991.86795711963, 0., 0.)
    end
end
NonlinearBinaryFluid(params) = NonlinearBinaryFluid(params.consvar, params)

const NBF   = NonlinearBinaryFluid
const NBFS  = NonlinearBinaryFluid{:entropy}
const NBFP  = NonlinearBinaryFluid{:potential_temperature}
# const NBFC  = NonlinearBinaryFluid{:conservative_temperature}

const PVCons = NamedTuple{(:p, :v, :consvar)}


@fastmath @muladd @inlineall begin
    ## private functions
    function solve_quad(a, b, c, prec)
        # rudimentary function assumes a != 0, b >= 0 and discriminant >= 0
        # selects root with smallest magnitude
        δ = eps(prec)

        q = -0.5 * (b + sqrt(b^2 - 4a*c))
        
        # choose root closest to zero
        r1 = q / (a+δ)
        r2 = c / q
        R1 = abs(r1)
        R2 = abs(r2)
        m = (sign(R1 - R2) + 1)/2
        soln = r1 * (1 - m) + r2 * m

        return soln
    end

    # T(p, θ)
    NBF_temperature_from_pθ((; v0, Cp, α_T, γ, α_TT, p0, T0), p, θ) = (θ / Cp) * NBF_exner_from_pθ((; Cp, v0, α_T, γ, p0, α_TT, T0), p, θ)
    # θ(p, T)
    function NBF_θ_from_pT(fluid::NonlinearBinaryFluid{CV, F}, p, T) where {CV, F}
        (; v0, Cp, α_T, γ, α_TT, p0, T0) = fluid
        a = -(v0/Cp) * (p - p0) * α_TT
        b = 1. + (v0/Cp) * ( α_T * (1 + 0.5 * γ * (p - p0)) + α_TT * T0 ) * (p - p0)
        c = (v0/Cp) * α_T * T0 * (1 + 0.5 * γ * (p - p0)) * (p - p0) - (T - T0)
        soln = solve_quad(a, b, c, F)
        return T0 + soln
    end
    # θ(s)
    NBF_θ_from_s((; T0, η0, Cp), s) = T0 * exp( (s - η0) / Cp)
    # s(θ)
    NBF_s_from_θ((; T0, η0, Cp), θ) = η0 + Cp * log(θ / T0)
    # θ(p, v, q)
    function NBF_θ_from_pvq(fluid::NonlinearBinaryFluid{CV, F}, p, v, q) where {CV, F}
        (; T0, p0, S0, v0, α_T, α_TT, α_S, α_p, γ) = fluid
        a = 0.5 * α_TT
        b = α_T * (1 + γ * (p - p0))
        c = - α_S * (q - S0) - α_p * (p - p0) - (v - v0)/v0
        soln = solve_quad(a, b, c, F)
        return T0 + soln
    end
    
    # p(v, θ, q)
    NBF_pressure_from_vθq((; p0, α_T, γ, T0, α_p, v0, α_S, S0, α_TT), v, θ, q) = p0 + inv( α_T * γ * (θ - T0) - α_p ) * ( (v - v0)/v0 + α_S * (q - S0) - α_T * (θ - T0) - 0.5 * α_TT * (θ - T0)^2 )
    # v(p, θ, q)
    NBF_specific_volume_from_pθq((; p0, T0, S0, v0, α_p, α_T, α_S, γ, α_TT), p, θ, q) = v0 * ( 1 + α_T * (1 + γ * (p - p0)) * (θ - T0) - α_S * (q - S0) - α_p * (p - p0) + 0.5 * α_TT * (θ - T0)^2 )
    # π = Cp.T / θ
    NBF_exner_from_pθ((; Cp, v0, α_T, γ, p0, α_TT, T0), p, θ) = Cp + v0 * ( α_T * (1 + 0.5 * γ * (p - p0)) + α_TT * (θ - T0) ) * (p - p0)
    
    function NBF_θ_from_vTq(fluid, v, T, q)
        (; v0, T0, S0, α_p, α_T, α_S, α_TT, γ, Cp) = fluid
        # constructed polynomial in Δθ by eliminating Δp in v(p, θ, q) and T(p, θ)
        
        # helper parameters
        A = (v - v0)/v0 + α_S * (q - S0)
        b = 0.5 * α_TT
        d = 0.5 * α_T * γ
        ΔT = T - T0

        # polynomial coefficients
        c0 = -Cp * α_p^2 * ΔT + v0 * (d * A^2 * T0 - α_T * α_p * A * T0)
        c1 = Cp * (α_p^2 + 4 * α_p * d * ΔT) + v0 * ( α_T^2 * α_p * T0 + d * A^2 - α_T * α_p * A - 2 * b * α_p * A * T0)
        c2 = -Cp * (4 * α_p * d + 4 * d^2 * ΔT) + v0 * (α_T^2 * α_p + 3 * α_T * b * α_p * T0 + 3 * b * d * T0 * A - α_T^2 * d * T0 - 2 * b * α_p * A)
        c3 = 4 * Cp * d^2 + v0 * (2 * b^2 * α_p * T0 + 2 * b * d * A - α_T^2 * d - 4 * α_T * b * d * T0 + 3 * α_T * b * α_p)
        c4 = v0 * (2 * b^2 * α_p - 4 * α_T * b * d - 3 * b^2 * d * T0)
        c5 = - 3 * b^2 * d * v0

        # solve polynomial P(θ; v, T, q) = 0
        rts = PolynomialRoots.roots([c0, c1, c2, c3, c4, c5])
        rts_real = real.(rts[iszero.(imag.(rts))])

        # discard complex roots and find solution closest to T0
        _, idx = findmin(abs.(rts_real))

        return T0 + rts_real[idx]
    end

    # h(p, θ, q)
    NBF_h_from_pθq((; p0, T0, S0, v0, α_p, α_T, α_S, γ, α_TT, Cp, η0, μ0), p, θ, q) = Cp * (θ - T0) + μ0 * (q - S0) + v0 * ( 1 + α_T * (1 + 0.5 * γ * (p - p0)) * (θ - T0) - α_S * (q - S0) - 0.5 * α_p * (p - p0) + 0.5 * α_TT * (θ - T0)^2 ) * (p - p0)
    # c2(p, θ, q)
    NBF_sound_speed2_from_pθq((; p0, T0, S0, v0, α_p, α_T, α_S, γ, α_TT), p, θ, q) = NBF_specific_volume_from_pθq((; p0, T0, S0, v0, α_p, α_T, α_S, γ, α_TT), p, θ, q)^2 / (v0 * (α_p - α_T * γ * (θ - T0)))
    # μ(p, T, q)
    NBF_chemical_potential_from_p((; μ0, v0, α_S, p0), p, q) = μ0 * one(q) - v0 * α_S * (p - p0)
    
    # S(p, v, θ)
    NBF_salinity_from_pvθ((; p0, T0, S0, v0, α_p, α_T, α_S, γ, α_TT), p, v, θ) = S0 + (-(v - v0)/v0 + α_T * (1 + γ * (p - p0)) * (θ - T0) - α_p * (p - p0) + 0.5 * α_TT * (θ - T0)^2 )/α_S

    ## all consvar
    # defined throughout ClimFluids:
    specific_entropy(           fluid::NBF, (p, T, q)::PTQ) = NBF_s_from_θ(fluid, NBF_θ_from_pT(fluid, p, T))
    specific_enthalpy(          fluid::NBF, (p, T, q)::PTQ) = NBF_h_from_pθq(fluid, p, NBF_θ_from_pT(fluid, p, T), q)
    specific_internal_energy(   fluid::NBF, (p, T, q)::PTQ) = NBF_h_from_pθq(fluid, p, NBF_θ_from_pT(fluid, p, T), q) - p * NBF_specific_volume_from_pθq(fluid, p, NBF_θ_from_pT(fluid, p, T), q)
    potential_temperature(      fluid::NBF, (p, T, q)::PTQ) = NBF_θ_from_pT(fluid, p, T)
    potential_enthalpy(         fluid::NBF, (p, T, q)::PTQ) = NBF_h_from_pθq(fluid, fluid.p0, NBF_θ_from_pT(fluid, p, T), q)
    potential_volume(           fluid::NBF, (p, T, q)::PTQ) = NBF_specific_volume_from_pθq(fluid, fluid.p0, NBF_θ_from_pT(fluid, p, T), q)
    heat_capacity(              fluid::NBF, (p, T, q)::PTQ) = fluid.Cp * one(T)
    pressure(                   fluid::NBF, (v, T, q)::VTQ) = NBF_pressure_from_vθq(fluid, v, NBF_θ_from_vTq(fluid, v, T, q), q)
    sound_speed2(               fluid::NBF, (p, T, q)::PTQ) = NBF_sound_speed2_from_pθq(fluid, p, NBF_θ_from_pT(fluid, p, T), q)
    
    # defined here only:
    potential_temperature(      fluid::NBF, (v, T, q)::VTQ) = NBF_θ_from_vTq(fluid, v, T, q)
    specific_volume(            fluid::NBF, (p, T, q)::PTQ) = NBF_specific_volume_from_pθq(fluid, p, NBF_θ_from_pT(fluid, p, T), q)
    chemical_potential(         fluid::NBF, (p, T, q)::PTQ) = NBF_chemical_potential_from_p(fluid, p, q)
    ds_dq(                      fluid::NBF, (p, T, q)::PTQ) = zero(T)   # ∂s/∂q with p, T const.
    
    # specific_entropy(           fluid::NBF, (p, θ, q)::PThQ) = NBF_s_from_θ(fluid, θ)

    ## consvar = :entropy
    
    # PTQ
    conservative_variable(  fluid::NBFS, (p, T, q)::PTQ)       = NBF_s_from_θ(fluid, NBF_θ_from_pT(fluid, p, T))
    conjugate_variable(     fluid::NBFS, (p, T, q)::PTQ)       = T
    modified_chemical_potential(fluid::NBFS, (p, T, q)::PTQ)   = NBF_chemical_potential_from_p(fluid, p, q)

    # PConsQ
    temperature(        fluid::NBFS, (p, s, q)::PConsQ)         = NBF_temperature_from_pθ(fluid, p, NBF_θ_from_s(fluid, s))
    specific_enthalpy(  fluid::NBFS, (p, s, q)::PConsQ)         = NBF_h_from_pθq(fluid, p, NBF_θ_from_s(fluid, s), q)
    specific_volume(    fluid::NBFS, (p, s, q)::PConsQ)         = NBF_specific_volume_from_pθq(fluid, p, NBF_θ_from_s(fluid, s), q)
    heat_capacity(      fluid::NBFS, (p, s, q)::PConsQ)         = fluid.Cp * one(s)
    modified_chemical_potential(fluid::NBFS, (p, s, q)::PConsQ) = NBF_chemical_potential_from_p(fluid, p, q)
    dcons_dq(               fluid::NBFS, (p, s, q)::PConsQ)     = zero(s) # with p, T const
    potential_temperature(  fluid::NBFS, (p, s, q)::PConsQ)     = NBF_θ_from_s(fluid, s)
    
    function exner_functions(fluid::NBFS, (p, s, q)::PConsQ)    # returns h, v, conjvar = T
        θ = NBF_θ_from_s(fluid, s)
        h = NBF_h_from_pθq(fluid, p, θ, q)
        v = NBF_specific_volume_from_pθq(fluid, p, θ, q)
        conjvar = NBF_temperature_from_pθ(fluid, p, θ)
        return h, v, conjvar
    end
    function volume_functions(fluid::NBFS, (p, s, q)::PConsQ)   # returns v, ∂v/∂p, ∂v/∂s, ∂v/∂q
        (; p0, v0, T0, α_p, α_T, α_S, γ, α_TT, Cp) = fluid
        θ = NBF_θ_from_s(fluid, s)
        v = NBF_specific_volume_from_pθq(fluid, p, θ, q)
        dv_dp = -v0 * α_p + v0 * α_T * γ * (θ - T0)
        dv_ds = (θ / Cp) * (v0 * α_T * ( 1 + γ * (p - p0)) + v0 * α_TT * (θ - T0))
        dv_dq = -v0 * α_S
        return v, dv_dp, dv_ds, dv_dq
    end

    # VCons
    temperature(    fluid::NBFS, (v, s, q)::VConsQ)    = NBF_temperature_from_pθ(fluid, NBF_pressure_from_vθq(fluid, v, NBF_θ_from_s(fluid, s), q), NBF_θ_from_s(fluid, s))
    pressure(       fluid::NBFS, (v, s, q)::VConsQ)    = NBF_pressure_from_vθq(fluid, v, NBF_θ_from_s(fluid, s), q)

    salinity(       fluid::NBFS, (p, v, s)::PVCons)      = NBF_salinity_from_pvθ(fluid, p, v, NBF_θ_from_s(fluid, s))


    ## consvar = :potential_temperature
    
    # PTQ
    conservative_variable(  fluid::NBFP, (p, T, q)::PTQ)       = NBF_θ_from_pT(fluid, p, T)
    conjugate_variable(     fluid::NBFP, (p, T, q)::PTQ)       = NBF_exner_from_pθ(fluid, p, NBF_θ_from_pT(fluid, p, T))
    modified_chemical_potential(fluid::NBFP, (p, T, q)::PTQ)   = NBF_chemical_potential_from_p(fluid, p, q)

    # PConsQ
    temperature(        fluid::NBFP, (p, θ, q)::PConsQ)        = NBF_temperature_from_pθ(fluid, p, θ)
    specific_enthalpy(  fluid::NBFP, (p, θ, q)::PConsQ)        = NBF_h_from_pθq(fluid, p, θ, q)
    specific_volume(    fluid::NBFP, (p, θ, q)::PConsQ)        = NBF_specific_volume_from_pθq(fluid, p, θ, q)
    heat_capacity(      fluid::NBFP, (p, θ, q)::PConsQ)        = fluid.Cp * one(θ)
    modified_chemical_potential(fluid::NBFP, (p, θ, q)::PConsQ)   = NBF_chemical_potential_from_p(fluid, p, q)
    dcons_dq(               fluid::NBFP, (p, θ, q)::PConsQ)   = zero(θ) # with p, T const
    potential_temperature(fluid::NBFP, (p, θ, q)::PConsQ)      = θ

    function exner_functions(fluid::NBFP, (p, θ, q)::PConsQ)    # returns h, v, conjvar = Cp * T / θ
        h = NBF_h_from_pθq(fluid, p, θ, q)
        v = NBF_specific_volume_from_pθq(fluid, p, θ, q)
        conjvar = fluid.Cp * NBF_temperature_from_pθ(fluid, p, θ) / θ
        return h, v, conjvar
    end
    function volume_functions(fluid::NBFP, (p, θ, q)::PConsQ)   # returns v, ∂v/∂p, ∂v/∂θ
        (; p0, v0, T0, α_p, α_T, α_S, γ, α_TT) = fluid
        v = NBF_specific_volume_from_pθq(fluid, p, θ, q)
        dv_dp = -v0 * α_p + v0 * α_T * γ * (θ - T0)
        dv_dθ = v0 * α_T * ( 1 + γ * (p - p0)) + v0 * α_TT * (θ - T0)
        dv_dq = -v0 * α_S
        return v, dv_dp, dv_dθ, dv_dq
    end

    # VCons
    temperature(    fluid::NBFP, (v, θ, q)::VConsQ)    = NBF_temperature_from_pθ(fluid, NBF_pressure_from_vθq(fluid, v, θ, q), θ)
    pressure(       fluid::NBFP, (v, θ, q)::VConsQ)    = NBF_pressure_from_vθq(fluid, v, θ, q)

    salinity(       fluid::NBFP, (p, v, θ)::PVCons)      = NBF_salinity_from_pvθ(fluid, p, v, θ)


    ## Fallback: convert to (p, T, q) if not implemented above
    canonical_state(fluid::NBF, (p, T, q)::PTQ)                 = (p = p, T = T, q = q)
    canonical_state(fluid::NBF, (p, consvar, q)::PConsQ)        = (p, T = NBF_temperature_from_pθ(fluid, p, potential_temperature(fluid, (; p, consvar, q))), q)
    canonical_state(fluid::NBF, (p, s, q)::PSQ)                 = (p, T = NBF_temperature_from_pθ(fluid, p, NBF_θ_from_s(fluid, s)), q)
    canonical_state(fluid::NBF, (p, θ, q)::PThQ)                = (p, T = NBF_temperature_from_pθ(fluid, p, θ), q)
    canonical_state(fluid::NBF, (p, v, q)::PVQ)                 = (p, T = NBF_temperature_from_pθ(fluid, p, NBF_θ_from_pvq(fluid, p, v, q)), q)
    canonical_state(fluid::NBF, (v, T, q)::VTQ)                 = (p = pressure(fluid, (; v, T, q)), T, q)
    canonical_state(fluid::NBF, (v, s, q)::VSQ)                 = canonical_state_NBF_vθ(fluid, (v, NBF_θ_from_s(fluid, s), q))
    canonical_state(fluid::NBFS, (v, s, q)::VConsQ)             = canonical_state_NBF_vθ(fluid, (v, NBF_θ_from_s(fluid, s), q))
    canonical_state(fluid::NBFP, (v, θ, q)::VConsQ)             = canonical_state_NBF_vθ(fluid, (v, θ, q))

    function canonical_state_NBF_vθ(fluid, (v, θ, q))
        p = NBF_pressure_from_vθq(fluid, v, θ, q)
        T = NBF_temperature_from_pθ(fluid, p, θ)
        return (; p, T, q)
    end

    exner_functions(fluid::NBF, (p, T, q)::PTQ)     = exner_functions(fluid,   (; p, consvar=conservative_variable(fluid, (; p, T, q)), q))
    volume_functions(fluid::NBF, (p, T, q)::PTQ)    = volume_functions(fluid,  (; p, consvar=conservative_variable(fluid, (; p, T, q)), q))
end