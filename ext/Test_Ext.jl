module Test_Ext

import ClimFluids: test_fluid
using ClimFluids
using Test

# extend isapprox to tuples, but only locally
≈(a, b) = isapprox(a, b)
≈(a::T, b::T) where {T<:Tuple} = all(map(isapprox, a, b))

function test_fluid(fluid, p, T)
    @testset "$(typeof(fluid))" begin
        test_state_functions(fluid, p, T)
        test_exner_functions(fluid, p, T)
        test_volume_functions(fluid, p, T)
    end
end

function test_exner_functions(fluid, p, T)
    fluid_pT = fluid(:p, :T)
    consvar = fluid_pT.conservative_variable(p, T)
    enthalpy(p, consvar) = ClimFluids.specific_enthalpy(fluid, (; p, consvar))
    h, v, exner = ClimFluids.exner_functions(fluid, (;p, consvar))
    h_, v_, exner_ = ClimFluids.fwdd_exner_functions(fluid, (;p, consvar))
    @test h ≈ enthalpy(p, consvar)
    @test v ≈ ClimFluids.specific_volume(fluid, (;p, consvar))
    @test h ≈ h_
    @test v ≈ v_
    @test exner ≈ exner_
end

function test_volume_functions(fluid, p, T)
    fluid_pT = fluid(:p, :T)
    consvar = fluid_pT.conservative_variable(p, T)
    volume(p, consvar) = ClimFluids.specific_volume(fluid, (; p, consvar))
    v, dv_dp, dv_dconsvar = ClimFluids.volume_functions(fluid, (;p, consvar))
    v_, dv_dp_, dv_dconsvar_ = ClimFluids.fwdd_volume_functions(fluid, (;p, consvar))
    @test v ≈ volume(p, consvar)
    @test dv_dp ≈ dv_dp_
    @test dv_dconsvar ≈ dv_dconsvar_
end

function test_state_functions(fluid, p, T)
    # check that all state functions can be computed
    # from all combinations of state variables (p or v, T or S or consvar, q)
    fluid_pT = fluid(:p, :T)
    v = fluid_pT.specific_volume(p, T)
    consvar = fluid_pT.conservative_variable(p, T)
    s = fluid_pT.specific_entropy(p, T)
    states = ((; p, T), (; p, s), (; p, consvar), (; v, T), (; v, s), (; v, consvar))
    for fun in ClimFluids.all_state_functions()
        @testset "$fun(::$(typeof(fluid)))" begin
            state = (; p, v)
            val = fun(fluid, state)
            # @btime $fun($fluid, $state)
            for state in states
                @test fun(fluid, state) ≈ val
                # @btime $fun($fluid, $state)
            end
        end
    end
end

end # module
