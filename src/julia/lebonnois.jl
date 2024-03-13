#===============================================================================#
# Ajout de d'un gaz avec Cp varible  (Work in progress)                         #
#-------------------------------------------------------------------------------#

struct VarCpPerfectGas{F} <: PerfectGas{F}
  kappa0::F
  Cp0   ::F
  Cv0   ::F
  R     ::F
  p0    ::F
  T0    ::F
  nu    ::F
  inv_T0 ::F
  inv_p0 ::F
  inv_Cp0 ::F
  inv_nu :: F
  inv_nup1 :: F
end
VCPG{F} = VarCpPerfectGas{F}

"""
    VarCpPerfectGas(kappa,Cp,p0,T0,nu)
Returns an object describing the thermodynamics of a single-component perfect gas with variable Cp. 
This object can then be used as first argument of thermodynamic functions.
"""
VarCpPerfectGas(kappa0, Cp0, p0, T0, nu) = VarCpPerfectGas(
    kappa0, Cp0, (1-kappa0)*Cp0, kappa0*Cp0, p0, T0, nu, 
    inv(T0), inv(p0), inv(Cp0), inv(nu), inv(1+nu))
VarCpPerfectGas(params) =  VarCpPerfectGas(
    params.kappa0, params.Cp0, params.p0, params.T0, params.nu)


@fastmath @muladd @inlineall begin
    # private functions
    LB10_enthalpy(Cp0,T0,nu, T) = inv(1+nu)*Cp0*T0*powm1(T*inv(T0), 1+nu)

    # interface functions
    function heat_capacity(gas::VCPG, (p,T)::PT) 
        (;nu, T0, p0, kappa0, Cp0, R)=gas
        return Cp0 * (T*inv(T0))^nu
    end

    specific_enthalpy(gas::VCPG, (p,T)::PT) = LB10_enthalpy(gas.Cp0,gas.T0,gas.nu, T) 

    # avoid loss of accuracy when nu is small (IPG limit)
    powm1(x,y) = expm1(y*log(x)) # x^y-1, accurate for small y
    pow1p(x,y) = exp(y*log1p(x)) # (1+x)^y, accurate for small x
    pow1pm1(x,y) = expm1(y*log1p(x)) # (1+x)^y-1, accurate for small x,y

    function temperature(gas::VCPG, (p,theta)::PTh) 
        (;nu, T0, p0, kappa0, Cp0, R)=gas
        # return (theta^nu + kappa0*nu*T0^nu*log(inv(p0)*p) )^inv(nu)
        # X = (theta/T0)^nu-1
        X = powm1(theta*inv(T0),nu)
        # Y = ((theta/T0)^nu+nu.kappa.log(p/p0))^(1/nu)
        Y = pow1p( X + nu*kappa0*log(p*inv(p0)), inv(nu) ) 
        return T0*Y 
    end
    function potential_temperature(gas::VCPG, (p,T)::PT)
        (;nu, T0, p0, kappa0, Cp0, R)=gas
        # return (T^nu - kappa0*nu*T0^nu*log(inv(p0)*p) )^inv(nu)
        # X = (T/T0)^nu-1
        X = powm1(T*inv(T0),nu) 
        # Y = ((T/T0)^nu-nu.kappa.log(p/p0))^(1/nu)
        Y = pow1p( X-nu*kappa0*log(inv(p0)*p), inv(nu) ) 
        return T0*Y 
    end
    
    conservative_variable(gas::VCPG, (p,T)::PT) = specific_entropy(gas, (;p,T))
    function specific_entropy(gas::VCPG, (p,T)::PT)
        (; nu, T0, p0, Cp0, R)=gas
        T,p = T*inv(T0), p*inv(p0)
        return (Cp0*inv(nu))*powm1(T,nu)- R*log(p)
    end
    function temperature(gas::VCPG, (p,s)::PCons)
        # s = (Cp0/nu)*((T/T0)^nu-1) - R*log(p/p0)
        # (T/T0)^nu = 1 + nu*(s/Cp0+kappa0*log(p/p0))
        # T = T0*(1+nu*(s/Cp0+kappa0*log(p/p0)))^(1/nu)
        (;nu, T0, p0, kappa0, Cp0)=gas
        X = s*inv(Cp0)+kappa0*log(p*inv(p0))
        return T0*pow1p(nu*X, inv(nu))
    end
    function exner_functions(gas::VCPG, (p,s)::PCons)
        (;nu, T0, p0, kappa0, Cp0, R)=gas
        X = s*inv(Cp0)+kappa0*log(p*inv(p0))
        Z = pow1pm1(nu*X, inv(nu)) # (1+nuX)^(1/nu)-1
        T = T0*(1+Z)   
        h = inv(1+nu)*Cp0*T0*( Z+nu*X*(1+Z) )
        return h-T*s, s, T
    end

    function exner_functions(gas::VCPG, (p,T)::PT)
        consvar = conservative_variable(gas, (;p,T))
        exner_functions(gas, (; p,consvar))
    end

    function exner_functions(gas::VCPG, (p,theta)::PTh)
        (;nu, T0, kappa0, Cp0, inv_p0, inv_Cp0, inv_nu, inv_nup1) = gas
        X = s*inv_Cp0+kappa0*log(p*inv_p0)
        Z = pow1pm1(nu*X, inv_nu) # (1+nuX)^(1/nu)-1
        T = T0*(1+Z)
        h = inv_nup1*Cp0*T0*( Z+nu*X*(1+Z) )
        return h-T*s, s, T
    end

    # canonical_state : p,T
    canonical_state(gas::VCPG, (p,v)::PV)          = (; p, T=p*v/gas.R)
    canonical_state(gas::VCPG, (p,consvar)::PCons) = (; p, T=temperature(gas, (;p,consvar)))
    canonical_state(gas::VCPG, (p,s)::PS)          = (; p, T=temperature(gas, (;p,consvar=s)))
    canonical_state(gas::VCPG, (v,T)::VT)    = (; p=gas.R*T/v, T)
    canonical_state(gas::VCPG, state::VS)    = canonical_state_VCPG_vs(gas, state)
    canonical_state(gas::VCPG, state::VCons) = canonical_state_VCPG_vs(gas, state)

    function canonical_state_VCPG_vs(gas, (v,s))
        (; R, T0, Cp0, inv_p0, nu) = gas
        # Newton-Raphson solving for T :
        #  0 = f(T) = (Cp0/nu)*((T/T0)^nu-1) - R*log(RT/p0.v0) + R*log(v/v0) - s
        #     df/dT = (Cp0/T)*exp(nu*log(T/T0)) - R/T
        # Initial guess as in an IPG (valid for small nu)
        p, T = canonical_state_IPG_vs((; T0, R, inv_p0, Cv=Cp0-R), (v,s))
        for _ =  1:5
            f = (Cp0/nu)*powm1(T/T0, nu) - R*log(R*T*inv_p0) + R*log(v) - s
            Tdf = Cp0*exp(nu*log(T/T0))-R
            T = T*(1-f/Tdf)
        end
        return (p=R*T/v, T)
    end

end # @fastmath @muladd
