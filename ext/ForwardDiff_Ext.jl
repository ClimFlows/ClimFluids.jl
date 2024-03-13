using .ForwardDiff: Partials, Dual

@inline duals(vars...) = duals_(vars, Val(length(vars)))
@inline duals_(vars, N) = ntuple(i->dual(vars[i], N,i), N)
@inline dual(var,N,i)   = Dual(var, Partials(ntuple(j->dirac(var,i,j), N)))
@inline dirac(var,i,j)  = (j==i) ? one(var) : zero(var)

@inline function exner_fwdd(gas::SimpleFluid, (p, consvar_)::PCons)
    consvar = Dual(consvar_, one(consvar_))
    T = temperature(gas, (;p, consvar))
    h = specific_enthalpy(gas, (;p,T))
    (hh, (Pi0,)) = h.value, h.partials
    return hh-Pi0*consvar_, consvar_, Pi0
end

@inline function exner_fwdd(gas::BinaryFluid, (p, consvar_, q_)::PConsQ)
    consvar, q = duals(consvar_, q_)
    T = temperature(gas, (;p, consvar, q))
    h = specific_enthalpy(gas, (;p,T,q))
    (hh, (Pi0, Pi1)) = h.value, h.partials
    return hh-Pi0*consvar_-Pi1*q_, consvar_, Pi0, q_, Pi1
end
