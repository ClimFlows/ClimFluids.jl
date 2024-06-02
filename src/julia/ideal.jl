"""
    fluid = IdealPerfectGas(consvar, kappa, Cp, p0, T0)
    fluid = IdealPerfectGas(options)
Return an object `fluid` describing the thermodynamics of a single-component ideal perfect gas.
`fluid` can then be used as first argument of thermodynamic functions.
Argument `consvar` defines the conservative variable chosen. Valid values
are `:temperature`, `:entropy`, `:enthalpy`, `:volume`, corresponding to
potential temperature, entropy, potential enthalpy, potential specific volume.

Arguments can also be passed as a named tuple `options`.
"""
struct IdealPerfectGas{CV,F} <: PerfectGas{F}
    kappa::F
    Cp::F
    Cv::F
    R ::F
    p0::F
    T0::F
    inv_p0::F
    inv_Cp::F

    function IdealPerfectGas(consvar::Symbol, kappa,Cp,p0,T0)
        @assert consvar in (:temperature, :entropy, :enthalpy, :volume)
        new{consvar,typeof(kappa)}(kappa, Cp, (1-kappa)*Cp, kappa*Cp, p0, T0, inv(p0), inv(Cp))
    end
end
IdealPerfectGas(params) = IdealPerfectGas(params.consvar,params.kappa, params.Cp, params.p0, params.T0)

const IPG = IdealPerfectGas
const IPGT = IdealPerfectGas{:temperature}
const IPGS = IdealPerfectGas{:entropy}
const IPGH = IdealPerfectGas{:enthalpy}
const IPGV = IdealPerfectGas{:volume}

@fastmath @muladd @inlineall begin
    # private routines
    IPG_enthalpy((;Cp), T)             = Cp*T
    IPG_entropy((;Cp,R,p0,T0), p,T)    = Cp*log(T*inv(T0))-R*log(p*inv(p0))
    IPG_temperature_ps((;inv_Cp,kappa,p0,T0), p,s)= T0*exp(s*inv_Cp+kappa*log(p*inv(p0)))
    IPG_internal_energy((;Cv), T)      = Cv*T
    IPG_sound_speed2((;R, kappa), T)   = R*inv(1-kappa)*T
    IPG_exner((;kappa, inv_p0), p)     = (p*inv_p0)^kappa
    IPG_theta((;kappa, inv_p0), p, T)  = T*(p*inv_p0)^-kappa
    IPG_pot_volume(gas, p,T) = (gas.R*gas.inv_p0)*IPG_theta(gas, p,T)
    # routines independent from choice of conservative variable
    specific_entropy(        gas::IPG, (p,T)::PT)     = IPG_entropy(gas, p,T)
    specific_enthalpy(       gas::IPG, (p,T)::PT)     = IPG_enthalpy(gas, T)
    specific_internal_energy(gas::IPG, (p,T)::PT)     = IPG_internal_energy(gas,T)
    sound_speed2(            gas::IPG, (p,T)::PT)     = IPG_sound_speed2(gas, T)
    potential_temperature(   gas::IPG, (p,T)::PT)     = IPG_theta(gas, p, T)
    potential_enthalpy(      gas::IPG, (p,T)::PT)     = IPG_enthalpy(gas, IPG_theta(gas, p, T))
    potential_volume(        gas::IPG, (p,T)::PT)     = IPG_pot_volume(gas, p, T)
    temperature(             gas::IPG, (p,s)::PS)     = IPG_temperature_ps(gas, p, s)
    exner_functions(         gas::IPG, (p,T)::PT)     = exner_functions(gas, (; p, consvar=conservative_variable(gas, (;p,T))))
    # routines depending on conservative variable
    #   potential temperature
    conservative_variable(   gas::IPGT, (p,T)::PT) = IPG_theta(gas, p,T)
    temperature(             gas::IPGT, (p,theta)::PCons) = theta*IPG_exner(gas, p)
    exner_functions(gas::IPGT, (p,theta)::PCons) = let Π = gas.Cp*IPG_exner(gas, p), h = theta*Π
        h, gas.kappa*h*inv(p), Π
    end
    #   potential enthalpy
    conservative_variable(   gas::IPGH, (p,T)::PT) = gas.Cp*IPG_theta(gas, p,T)
    temperature(             gas::IPGH, (p,hpot)::PCons) = gas.inv_Cp*hpot*IPG_exner(gas, p)
    exner_functions(         gas::IPGH, (p,hpot)::PCons) = let exner = IPG_exner(gas, p), h = hpot*exner
        h, gas.kappa*h*inv(p), exner
    end
    #   entropy
    conservative_variable(gas::IPGS, (p,T)::PT) = IPG_entropy(gas, p,T)
    function temperature(gas::IPGS, (p,s)::PCons)
        (; inv_Cp, T0, R, inv_p0) = gas
        return T0*exp(inv_Cp*(s+R*log(p*inv_p0)))
    end
    function exner_functions(gas::IPGS, (p,s)::PCons)
        T = temperature(gas, (;p,consvar=s))
        return gas.Cp*T, gas.R*T*inv(p), T
    end

    # Allows fallback implementations for inputs not covered above. Converts to (p,T)
    canonical_state(gas::IPG, (p,consvar)::PCons) = (p, T=temperature(gas, (; p,consvar)))
    canonical_state(gas::IPG, (p,s)::PS) = (p, T=temperature(gas, (; p,s)))
    canonical_state(gas::IPG, (p,v)::PV) = (p, T=p*v/gas.R)
    canonical_state(gas::IPG, (v,T)::VT) = (p=gas.R*T/v, T)
    canonical_state(gas::IPG, (v,s)::VS) = canonical_state_IPG_vs(gas, (v,s))
    canonical_state(gas::IPGS, (v,s)::VCons) = canonical_state_IPG_vs(gas, (v,s))
    canonical_state(gas::IPGT, (v,theta)::VCons) = canonical_state_IPG_vtheta(gas, (v,theta))
    canonical_state(gas::IPGH, (v,hpot)::VCons) = canonical_state_IPG_vtheta(gas, (v,hpot/gas.Cp))
    function canonical_state_IPG_vs(gas, (v,s))
        # entropy = Cv*log(T*inv(T0))+R*log(v*inv(v0)) where v0=RT0/p0
        (; Cv, T0, R, inv_p0) = gas
        v0 = R*T0*inv_p0
        T = T0*exp((s-R*log(v/v0))/Cv)
        p = R*T/v
        return (p=R*T/v, T)
    end
    function canonical_state_IPG_vtheta(gas, (v,theta))
        # theta = T*(p/p0)^-(R/Cp)
        # pv = RT
        # R*theta = pv*(p/p0)^-(R/Cp)
        # X = R*theta/p0*(v^-R/Cp) = (pv/p0)^(Cp-R)/Cp
        # pv = p0 * X^(Cp/Cv)
        (; Cp, Cv, R, p0, inv_p0) = gas
        X = (R*theta*inv_p0)*v^(-R/Cp)
        pv = p0*X^(Cp/Cv)
        return (p=pv/v, T=pv/R)
    end

end # @fastmath @muladd @inlineall
