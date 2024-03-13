struct IdealBinaryGas{F} <: BinaryFluid{F}
    Cpd::F 
    Rd ::F
    Cpv::F 
    Rv ::F
    p0::F
    T0::F
    inv_p0::F
    inv_T0::F

    IdealBinaryGas(Cpd, Rd, Cpv, Rv, p0, T0) = new{typeof(Cpd)}(
        Cpd, Rd, Cpv, Rv, p0, T0, inv(p0), inv(T0))
        
end

function IdealBinaryGas(options)
    (;Cpd, Rd, Cpv, Rv, p0, T0) = options
    return IdealBinaryGas(Cpd, Rd, Cpv, Rv, p0, T0)
end

@fastmath @muladd @inlineall begin

    #== helper functions ==#

    mix(q,Xd,Xv) = (1-q)*Xd + q*Xv

    function kappa_Cp(q,Rd,Rv,Cpd,Cpv)
        Cp = mix(q,Cpd,Cpv)
        inv_Cp = inv(Cp)
        return mix(q,Rd,Rv)*inv(Cp), Cp, inv_Cp
    end

    function molar_ratio(Rd, Rv, q)
        # rNv = q / Rv
        # rNd = (1-q) / Rd
        # Xd = Nd/(Nv+Nd) = rNd/(rNv+rNd)
        #    = (1-q)*Rv / (q*Rd + (1-q)*Rv)
        # Xv = Nv/(Nv+Nd) = rNv/(rNv+rNd)
        #    = q*Rd / (q*Rd + (1-q)*Rv)
        X = inv(mix(q,Rv,Rd)) # this "inversion" of Rd and Rv is *correct*
        return X*Rv, X*Rd
        # NB : let Yd=X*Rv, Yv=X*Rd
        #      the molar ratios Nd/(Nd+Nv) and Nv/(Nd+Nv)
        #      are Xd=(1-q)*Yd and Xv=q*Yv
    end

    xlogx(x, C) = x*log(C*(x+eps(x))) # regularize x*log(C*x) for small, non-dimensional x

    #== thermodynamic functions ==#
    function specific_volume(gas::IdealBinaryGas, (p,T,q)::PTQ) 
        (; Rv, Rd) = gas
        return mix(q,Rd,Rv)*T*inv(p)
    end

    function temperature(gas::IdealBinaryGas, (p,v,q)::PVQ) 
        (; Rv, Rd) = gas
        return p*v*inv(mix(q,Rd,Rv))
    end

    function specific_enthalpy(gas::IdealBinaryGas, (p,T,q)::PTQ) 
        (; Cpd, Cpv, T0) = gas
        return mix(q,Cpd,Cpv)*T
    end


    function specific_entropy(gas::IdealBinaryGas, (p,T,q)::PTQ) 
        (; Rd, Rv, Cpd, Cpv, inv_T0, inv_p0) = gas
        Cp = mix(q,Cpd,Cpv)
        Xd, Xv = molar_ratio(Rd, Rv, q)
        return Cp*log(T*inv_T0)-Rv*xlogx(q,Xv*p*inv_p0)-Rd*xlogx(1-q, Xd*p*inv_p0)
    end

    conservative_variable(gas::IdealBinaryGas, pTq::PTQ) = potential_enthalpy(gas, pTq)

    function potential_temperature(gas::IdealBinaryGas, (p,T,q)::PTQ) 
        (; Rd, Rv, Cpd, Cpv, inv_T0, inv_p0) = gas
        kappa = mix(q,Rd,Rv)*inv(mix(q,Cpd,Cpv))
        return T*(p*inv_p0)^(-kappa)
    end

    exner_factor(kappa, inv_p0, p) = pow_fast(p*inv_p0, kappa)

    function potential_enthalpy(gas::IdealBinaryGas, (p,T,q)::PTQ) 
        (; Rd, Rv, Cpd, Cpv, inv_T0, inv_p0) = gas
        kappa, Cp, _ = kappa_Cp(q,Rd,Rv,Cpd,Cpv)
        theta = T*exner_factor(-kappa, inv_p0, p)
        return Cp*theta
    end

    function temperature(gas::IdealBinaryGas, (p,hpot,q)::PConsQ) 
        (; Rd, Rv, Cpd, Cpv, T0, inv_p0) = gas
        kappa, _, inv_Cp = kappa_Cp(q,Rd,Rv,Cpd,Cpv)
        theta = inv_Cp*hpot
        return theta*exner_factor(kappa, inv_p0, p)
    end

    function exner_functions(gas::IdealBinaryGas, (p,hpot,q)::PConsQ) 
        (; Rd, Rv, Cpd, Cpv, T0, inv_p0) = gas
        kappa, Cp, inv_Cp = kappa_Cp(q,Rd,Rv,Cpd,Cpv)
        log_pi0 = kappa*log(p*inv_p0)
        pi0 = exp(log_pi0)
        T = inv_Cp*hpot*pi0
        h = Cp*T
        pi1 = log_pi0*T * ((Rv-Rd)*inv(kappa)+Cpd-Cpv)
        return h-pi0*hpot-pi1*q, hpot, pi0, q, pi1
    end

end # @fastmath @muladd @inlineall