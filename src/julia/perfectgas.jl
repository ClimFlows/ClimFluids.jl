abstract type PerfectGas{F} <: AbstractFluid{F,1} end

@fastmath @muladd @inlineall begin

    temperature(gas::PerfectGas, (p,v)::PV) = p*v*inv(gas.R)
    specific_volume(gas::PerfectGas, (p,T)::PT) = gas.R*T*inv(p)
    specific_internal_energy(gas::PerfectGas, (p,T)::PT) = specific_enthalpy(gas,p,T) - gas.R*T
    function sound_speed2(gas::PerfectGas, (p,T)::PT) 
        (R, Cp) = (gas.R, heat_capacity(gas, (;p,T)))
        return (Cp * inv(Cp-R) * R) * T
    end

end # @fastmath
