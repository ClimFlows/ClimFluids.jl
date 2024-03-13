import ClimFluids as Thermo
using Test
# using BenchmarkTools

Base.isapprox(a::T, b::T) where {T<:Tuple} = all(map(isapprox, a,b))

function test_fluid(fluid, p, T)
    testname = "$(typeof(fluid))"
    @testset "$testname" begin
        fluid_pT = fluid(:p, :T)
        v = fluid_pT.specific_volume(p,T)
        consvar = fluid_pT.conservative_variable(p,T)
        s = fluid_pT.specific_entropy(p,T)
        states = ((; p, T), (; p, s), (; p, consvar),
                (; v, T), (; v, s), (; v, consvar))
        for fun in Thermo.all_state_functions()
            @testset let name="$fun(::$(typeof(fluid)))"
                state = (; p,v)
                val = fun(fluid, state)
#                @btime $fun($fluid, $state)
                for state in states
                    @test fun(fluid, state) ≈ val
#                    @btime $fun($fluid, $state)
                end
            end
        end
    end
end

function test()
    params = (kappa=2/7, Cp=1000, p0=1e5, T0=300, nu=0.1)
    params = (kappa0=params.kappa, Cp0=params.Cp, params...)

    p,T = 1.1e5, 275
    IPG = Thermo.IdealPerfectGas
    CPV = Thermo.VarCpPerfectGas
    consvars = Dict( IPG=>(:temperature,:entropy), CPV=>(:temperature,) )

    params = map(Float32, params)
    p, T = map(Float32, (p,T))

    for Fluid in keys(consvars)
        for consvar in consvars[Fluid]
            test_fluid(Fluid((; consvar, params...)), p, T)
        end
    end
end

test()

