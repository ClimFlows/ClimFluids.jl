
# @inline exner_functions(gas::AbstractFluid, args...) = exner_zyg(gas, args...)

@inline function exner_zyg(gas::SimpleFluid, (p,consvar)::PCons)
    h, (Pi,) = Zygote.withgradient(consvar) do consvar
        T = temperature(gas, (;p, consvar))
        specific_enthalpy(gas, (;p,T))
    end
    return h-Pi*consvar, consvar, Pi
end

@inline function exner_zyg(gas::BinaryFluid, (p, consvar, q)::PConsQ)
    h, (Pi0, Pi1) = Zygote.withgradient(consvar,q) do consvar,q
        T = temperature(gas, (;p, consvar, q))
        specific_enthalpy(gas, (;p,T,q))
    end
    return h-Pi0*consvar-Pi1*q, consvar, Pi0, q, Pi1
end
