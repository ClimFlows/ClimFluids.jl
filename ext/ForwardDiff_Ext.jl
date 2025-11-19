module ForwardDiff_Ext

using ForwardDiff: Partials, Dual
using ClimFluids: PCons, PConsQ, SimpleFluid, BinaryFluid, temperature
using ClimFluids: specific_enthalpy, specific_volume
import ClimFluids: fwdd_exner_functions, fwdd_volume_functions

@inline duals(vars...) = duals_(vars, Val(length(vars)))
@inline duals_(vars, N) = ntuple(i -> dual(vars[i], N, i), N)
@inline dual(var, N, i) = Dual(var, Partials(ntuple(j -> dirac(var, i, j), N)))
@inline dirac(var, i, j) = (j == i) ? one(var) : zero(var)

@inline function fwdd_exner_functions(gas::SimpleFluid, (p_, consvar_)::PCons)
    p, consvar = duals(p_, consvar_)
    T = temperature(gas, (; p, consvar))
    h = specific_enthalpy(gas, (; p, T))
    (hh, (v, Pi0,)) = h.value, h.partials
    return hh, v, Pi0
end

@inline function fwdd_exner_functions(gas::BinaryFluid, (p_, consvar_, q_)::PConsQ)
    p, consvar, q = duals(p_, consvar_, q_)
    T = temperature(gas, (; p, consvar, q))
    h = specific_enthalpy(gas, (; p, T, q))
    (hh, (v, Pi0, Pi1)) = h.value, h.partials
    return hh, v, Pi0, Pi1
end

@inline function fwdd_volume_functions(gas::SimpleFluid, (p_, consvar_)::PCons)
    p, consvar = duals(p_, consvar_)
    T = temperature(gas, (; p, consvar))
    v = specific_volume(gas, (; p, T))
    (v, (dv_dp, dv_dconsvar,)) = v.value, v.partials
    return v, dv_dp, dv_dconsvar
end

@inline function fwdd_volume_functions(gas::BinaryFluid, (p_, consvar_, q_)::PConsQ)
    p, consvar, q = duals(p_, consvar_, q_)
    T = temperature(gas, (; p, consvar, q))
    v = specific_volume(gas, (; p, T, q))
    (v, (dv_dp, dv_dconsvar, dv_dq)) = v.value, v.partials
    return v, dv_dp, dv_dconsvar, dv_dq
end

end # module
