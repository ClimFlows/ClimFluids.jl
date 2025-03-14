using Test
using ForwardDiff
using ClimFluids

function test()
    params = (kappa = 2 / 7, Cp = 1000, p0 = 1e5, T0 = 300, nu = 0.1)
    params = (kappa0 = params.kappa, Cp0 = params.Cp, params...)

    p, T = 1.1e5, 275
    IPG = ClimFluids.IdealPerfectGas
    CPV = ClimFluids.VarCpPerfectGas
    LSF = ClimFluids.LinearSimpleFluid
    consvars = Dict(IPG => (:temperature, :entropy, :enthalpy), CPV => (:temperature,), LSF => (:entropy, :conservative_temperature))

    prec = Float64
    params = map(prec, params)
    params = (; prec, params...) # LSF
    p, T = map(prec, (p, T))

    for Fluid in keys(consvars)
        for consvar in consvars[Fluid]
            ClimFluids.test_fluid(Fluid((; consvar, params...)), p, T)
        end
    end
end


test()

