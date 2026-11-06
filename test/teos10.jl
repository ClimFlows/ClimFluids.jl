using Test, ForwardDiff, GibbsSeaWater
using ClimFluids
using ClimFluids: TEOS10Seawater

# GSW works in absolute salinity S_A [g/kg], Celsius temperatures [°C] and sea pressure [dbar]
const T0 = 273.15
pressure_Pa(p_dbar) = 101325 + 1e4 * p_dbar
derivative(f, x) = ForwardDiff.derivative(f, x)
allocated(fun, fluid, state) = (fun(fluid, state); @allocated fun(fluid, state))

# a grid over the range of validity of the 75-term polynomial: S_A [g/kg], Θ [°C], p - p0 [dbar]
const GRID = Iterators.product(range(0.5, 42, 7), range(-2, 40, 7), range(0, 6000, 5))
# the whole grid in Float64, its corners and centre in other float types, whose arithmetic is slower
grid(::Type{Float64}) = GRID
grid(::Type) = Iterators.product((0.5, 21.0, 42.0), (-2.0, 19.0, 40.0), (0.0, 3000.0, 6000.0))

const CONSVARS = (:conservative_temperature, :potential_temperature, :entropy)

@testset "TEOS10Seawater" begin

    @testset "check values" begin
        fluid = TEOS10Seawater(:conservative_temperature)
        # IAPWS SR7-09, table 6: pure water at p0, T = 273.15 K and 313.15 K, to the 9 digits given
        for (T, g, g_T, g_TT) in ((273.15, 0.101342743e3, 0.147644587, -0.154472324e2),
                                  (313.15, -0.116198898e5, -0.572365181e3, -0.133463968e2))
            gibbs = ClimFluids.TEOS_gibbs0(fluid, T, 0.0)
            @test gibbs[1] ≈ g rtol = 1e-8
            @test gibbs[2] ≈ g_T rtol = 1e-8
            @test gibbs[4] ≈ g_TT rtol = 1e-8
        end
        # Roquet et al. (2015), appendix C.3: S_A = 30 g/kg, Θ = 10 °C, p - p0 = 1000 dbar
        @test fluid(:p, :consvar, :q).specific_volume(pressure_Pa(1000), T0 + 10, 0.030) ≈ 9.732819628e-4 rtol = 1e-9
    end

    # The polynomials match GSW to round-off: errors within a few eps times the scale of each quantity, apart from the
    # dynamic enthalpy, which GSW evaluates from its own coefficients, rounded to 12 digits (ours are those of v divided
    # by k + 1), so that h, ∂h/∂Θ and ∂h/∂q differ from GSW by up to 10⁻¹² of v (p - p0). Derivatives at fixed φ are
    # brought to fixed Θ with ∂Θ/∂φ and ∂Θ/∂q at fixed φ.
    @testset "against GSW: $consvar" for consvar in CONSVARS
        fluid = TEOS10Seawater(consvar)
        (; Cp0, p0, q_star) = fluid
        fluid_pTq, fluid_pφq = fluid(:p, :T, :q), fluid(:p, :consvar, :q)
        for (S, Θ, p_dbar) in GRID
            p, q = pressure_Pa(p_dbar), S / 1000
            θ = gsw_pt_from_ct(S, Θ)
            φ = fluid(:p, :theta, :q).conservative_variable(p, T0 + θ, q)
            Θ_φ = derivative(φ -> fluid_pφq.potential_enthalpy(p, φ, q), φ) / Cp0
            Θ_q = derivative(q -> fluid_pφq.potential_enthalpy(p, φ, q), q) / Cp0
            v, v_p, v_φ, v_q = fluid_pφq.volume_functions(p, φ, q)
            conjvar = fluid_pφq.conjugate_variable(p, φ, q)
            v_SA, v_CT, v_P = gsw_specvol_first_derivatives(S, Θ, p_dbar)
            h_SA, h_CT = gsw_enthalpy_first_derivatives(S, Θ, p_dbar)
            dynamic = v * (p - p0)                                # scale of the dynamic enthalpy
            # the 75-term polynomials
            @test v ≈ gsw_specvol(S, Θ, p_dbar) rtol = 8eps()
            @test v_p ≈ v_P rtol = 8eps()
            @test v_φ / Θ_φ ≈ v_CT atol = 8eps() * v / 40
            @test v_q - v_φ / Θ_φ * Θ_q ≈ 1000v_SA atol = 8eps() * v / q_star
            @test fluid_pφq.sound_speed2(p, φ, q) ≈ gsw_sound_speed(S, Θ, p_dbar)^2 rtol = 16eps()
            @test fluid_pφq.specific_enthalpy(p, φ, q) ≈ gsw_enthalpy(S, Θ, p_dbar) atol = 8eps() * Cp0 * 40 + 1e-12 * dynamic
            @test conjvar / Θ_φ ≈ h_CT atol = 8eps() * Cp0 + 1e-12 * dynamic / 40
            # (∂h/∂q)_Θ = (∂h/∂q)_φ - (∂h/∂Θ) (∂Θ/∂q)_φ, where (∂Θ/∂q)_φ comes from Gibbs terms of order Cp0 T0 unless φ = Θ
            conversion = consvar == :conservative_temperature ? 0.0 : 8eps() * Cp0 * T0
            @test fluid_pφq.modified_chemical_potential(p, φ, q) - conjvar / Θ_φ * Θ_q ≈ 1000h_SA atol = conversion + 1e-12 * dynamic / q_star
            # the Gibbs function at p0
            @test fluid_pφq.potential_temperature(p, φ, q) ≈ T0 + θ atol = 4eps() * T0
            @test fluid_pφq.specific_entropy(p, φ, q) ≈ gsw_entropy_from_ct(S, Θ) atol = 8eps() * Cp0
            # the temperature derives from h, so it differs from the full TEOS-10 by the error of the 75-term fit
            t = gsw_t_from_ct(S, Θ, p_dbar)
            @test fluid_pφq.temperature(p, φ, q) ≈ T0 + t atol = 1e-2
            p_dbar <= 100 && @test fluid_pφq.temperature(p, φ, q) ≈ T0 + t atol = 1e-4
            @test fluid_pφq.heat_capacity(p, φ, q) ≈ gsw_cp_t_exact(S, t, p_dbar) rtol = 2e-3
            @test fluid_pTq.chemical_potential_derivatives(p, T0 + t, q)[1] ≈ 1000 * gsw_gibbs(1, 0, 0, S, t, p_dbar) atol = 100
            @test fluid_pTq.chemical_potential(p, T0 + t, q) ≈ 1000 * gsw_gibbs(1, 0, 0, S, t, p_dbar) atol = 100
        end
    end

    # all state functions derive from h(p, η, q), so the identities of thermodynamics hold to round-off, in any float
    # type: the tolerances are those of Float64 scaled by eps(F) / eps(Float64), on part of the grid for BigFloat
    @testset "identities in $F: $consvar" for F in (Float64, BigFloat), consvar in CONSVARS
        fluid = TEOS10Seawater(consvar; F)
        fluid_pTq, fluid_pφq, fluid_pηq = fluid(:p, :T, :q), fluid(:p, :consvar, :q), fluid(:p, :s, :q)
        scale = eps(F) / eps(Float64)
        for (S, Θ, p_dbar) in grid(F)
            p, q = F(pressure_Pa(p_dbar)), F(S) / 1000
            φ = fluid(:p, :theta, :q).conservative_variable(p, fluid.T0 + F(gsw_pt_from_ct(S, Θ)), q)
            T = fluid_pφq.temperature(p, φ, q)
            η = fluid_pφq.specific_entropy(p, φ, q)
            h, v, conjvar = fluid_pφq.exner_functions(p, φ, q)
            @test v ≈ derivative(p -> fluid_pφq.specific_enthalpy(p, φ, q), p) rtol = 1e-10scale
            @test conjvar ≈ derivative(φ -> fluid_pφq.specific_enthalpy(p, φ, q), φ) rtol = 1e-10scale
            @test fluid_pφq.modified_chemical_potential(p, φ, q) ≈ derivative(q -> fluid_pφq.specific_enthalpy(p, φ, q), q) rtol = 1e-9scale
            @test T ≈ derivative(η -> fluid_pηq.specific_enthalpy(p, η, q), η) rtol = 1e-12scale
            @test fluid_pφq.potential_enthalpy(p, φ, q) ≈ fluid_pφq.specific_enthalpy(fluid.p0, φ, q) rtol = 1e-12scale atol = 1e-6scale
            @test fluid_pφq.specific_internal_energy(p, φ, q) ≈ h - p * v rtol = 1e-12scale
            @test fluid_pφq.sound_speed2(p, φ, q) ≈ -v^2 / derivative(p -> fluid_pηq.specific_volume(p, η, q), p) rtol = 1e-10scale
            @test fluid_pφq.heat_capacity(p, φ, q) ≈ T * derivative(T -> fluid_pTq.specific_entropy(p, T, q), T) rtol = 1e-9scale
            @test fluid_pφq.dcons_dq(p, φ, q) ≈ derivative(q -> fluid_pTq.conservative_variable(p, T, q), q) rtol = 1e-8scale atol = 1e-12scale
            # μ = ∂g/∂q at fixed (p, T), g = h - T η, and its derivatives
            μ, μ_p, μ_T, μ_q = fluid_pTq.chemical_potential_derivatives(p, T, q)
            gibbs(p, T, q) = fluid_pTq.specific_enthalpy(p, T, q) - T * fluid_pTq.specific_entropy(p, T, q)
            @test μ ≈ derivative(q -> gibbs(p, T, q), q) rtol = 1e-8scale
            @test μ_p ≈ derivative(q -> fluid_pTq.specific_volume(p, T, q), q) rtol = 1e-8scale
            @test μ_T ≈ -derivative(q -> fluid_pTq.specific_entropy(p, T, q), q) rtol = 1e-8scale
            @test μ_q ≈ derivative(q -> fluid_pTq.chemical_potential_derivatives(p, T, q)[1], q) rtol = 1e-8scale
            @test all(fluid_pφq.chemical_potential_derivatives(p, φ, q) .≈ (μ, μ_p, μ_T, μ_q))
            # the chemical potential alone, from (p, T, q) and (p, φ, q), as the first of its derivatives: μ is a sum of
            # Gibbs terms of order Cp0 T0, which may nearly cancel, and the compiler may round it differently in the two
            # functions
            @test fluid_pTq.chemical_potential(p, T, q) ≈ μ atol = 8eps(F) * fluid.Cp0 * fluid.T0
            @test fluid_pφq.chemical_potential(p, φ, q) ≈ fluid_pφq.chemical_potential_derivatives(p, φ, q)[1] atol = 8eps(F) * fluid.Cp0 * fluid.T0
        end
    end

    # consistency between the state variables, and with ForwardDiff
    @testset "ClimFluids.test_fluid: $consvar" for consvar in CONSVARS
        ClimFluids.test_fluid(TEOS10Seawater(consvar), pressure_Pa(1000), 285.0, 0.035)
        ClimFluids.test_fluid(TEOS10Seawater(consvar), pressure_Pa(5000), 275.0, 0.034)
    end

    # the Newton iterations converge to the round-off of the float type
    @testset "conversions between conservative variables in $F" for F in (Float64, BigFloat)
        fluid = TEOS10Seawater(:conservative_temperature; F)
        for (S, Θ, _) in grid(F), q in (F(S) / 1000, zero(F))
            θ = ClimFluids.TEOS_θ_Θq(fluid, fluid.T0 + F(Θ), q)
            @test ClimFluids.TEOS_Θ_θq(fluid, θ, q) ≈ fluid.T0 + F(Θ) rtol = 2eps(F)
            @test ClimFluids.TEOS_θ_ηq(fluid, ClimFluids.TEOS_η_θq(fluid, θ, q), q) ≈ θ rtol = 2eps(F)
        end
    end

    @testset "inversions in $F" for F in (Float64, BigFloat)
        fluid = TEOS10Seawater(:conservative_temperature; F)
        for (S, Θ, p_dbar) in grid(F)
            p, q, φ = F(pressure_Pa(p_dbar)), F(S) / 1000, fluid.T0 + F(Θ)
            v = fluid(:p, :consvar, :q).specific_volume(p, φ, q)
            # the error of p is the round-off of v amplified by 1 / |∂v/∂p| ≈ 2 × 10¹² Pa kg/m³
            @test fluid(:v, :consvar, :q).pressure(v, φ, q) ≈ p atol = 1e4 * eps(v) / 4e-13
            T = fluid(:p, :consvar, :q).temperature(p, φ, q)
            @test fluid(:p, :T, :q).conservative_variable(p, T, q) ≈ φ rtol = 4eps(F)
            state = ClimFluids.canonical_state(fluid, (; v, T, q))
            @test state.p ≈ p atol = 1e4 * eps(v) / 4e-13
            @test state.consvar ≈ φ rtol = 4eps(F)
        end
    end

    # (p, v, q) has two roots on either side of the density maximum: the warm one is returned
    @testset "inversion of v: the warm root" begin
        fluid = TEOS10Seawater(:conservative_temperature)
        p, q = pressure_Pa(0), 0.0
        v_cold = fluid(:p, :consvar, :q).specific_volume(p, T0 + 1, q)
        Θ = fluid(:p, :v, :q).conservative_variable(p, v_cold, q)
        @test Θ > T0 + 4
        @test fluid(:p, :consvar, :q).specific_volume(p, Θ, q) ≈ v_cold rtol = 4eps()
    end

    @testset "fresh water: $consvar" for consvar in CONSVARS
        fluid = TEOS10Seawater(consvar)
        p, θ = pressure_Pa(100), T0 + 10
        for q in (0.0, 1e-12)
            φ = fluid(:p, :theta, :q).conservative_variable(p, θ, q)
            fluid_pφq = fluid(:p, :consvar, :q)
            @test all(isfinite, (fluid_pφq.volume_functions(p, φ, q)..., fluid_pφq.exner_functions(p, φ, q)...,
                                 fluid_pφq.temperature(p, φ, q), fluid_pφq.heat_capacity(p, φ, q)))
        end
        # the chemical potential diverges as ln q and ∂μ/∂q as 1/q, so that the barodiffusion -μ_p/μ_q vanishes as q
        T = fluid(:p, :theta, :q).temperature(p, θ, 1e-9)
        _, μ_p₁, _, μ_q₁ = fluid(:p, :T, :q).chemical_potential_derivatives(p, T, 1e-9)
        _, μ_p₂, _, μ_q₂ = fluid(:p, :T, :q).chemical_potential_derivatives(p, T, 1e-10)
        @test (μ_p₁ / μ_q₁) / (μ_p₂ / μ_q₂) ≈ 10 rtol = 1e-3
    end

    @testset "beyond the range of the fit: $consvar" for consvar in CONSVARS
        # values stay finite and consistent, but are extrapolated
        fluid = TEOS10Seawater(consvar)
        p, T, q = pressure_Pa(12000), T0 + 45, 0.05
        φ = fluid(:p, :T, :q).conservative_variable(p, T, q)
        @test fluid(:p, :consvar, :q).temperature(p, φ, q) ≈ T rtol = 4eps()
        @test all(isfinite, fluid(:p, :consvar, :q).volume_functions(p, φ, q))
    end

    @testset "float type $F: $consvar" for F in (Float32, BigFloat), consvar in CONSVARS
        fluid, fluid64 = TEOS10Seawater(consvar; F), TEOS10Seawater(consvar)
        p, T, q = F(pressure_Pa(1000)), F(285), F(0.035)
        φ = fluid(:p, :T, :q).conservative_variable(p, T, q)
        p64, T64, q64 = Float64.((p, T, q))
        φ64 = fluid64(:p, :T, :q).conservative_variable(p64, T64, q64)
        values = (fluid(:p, :consvar, :q).volume_functions(p, φ, q)..., fluid(:p, :consvar, :q).exner_functions(p, φ, q)...,
                  fluid(:p, :consvar, :q).temperature(p, φ, q), fluid(:v, :consvar, :q).pressure(F(9.7e-4), φ, q),
                  fluid(:p, :T, :q).chemical_potential_derivatives(p, T, q)...)
        values64 = (fluid64(:p, :consvar, :q).volume_functions(p64, φ64, q64)..., fluid64(:p, :consvar, :q).exner_functions(p64, φ64, q64)...,
                    fluid64(:p, :consvar, :q).temperature(p64, φ64, q64), fluid64(:v, :consvar, :q).pressure(9.7e-4, φ64, q64),
                    fluid64(:p, :T, :q).chemical_potential_derivatives(p64, T64, q64)...)
        @test all(x -> x isa F, (φ, values...))
        # Float32 to its precision, BigFloat to that of Float64: the same decimal coefficients, rounded to each type
        @test all(isapprox.(values, values64; rtol = F == Float32 ? 1e-4 : 1e-12))
    end

    @testset "no allocations, inferred: $consvar" for consvar in CONSVARS
        fluid = TEOS10Seawater(consvar)
        p, T, q = pressure_Pa(1000), 285.0, 0.035
        φ = fluid(:p, :T, :q).conservative_variable(p, T, q)
        state, state_T, state_v = (; p, consvar = φ, q), (; p, T, q), (; v = 9.7e-4, consvar = φ, q)
        for (fun, s) in ((ClimFluids.volume_functions, state), (ClimFluids.exner_functions, state),
                         (ClimFluids.temperature, state), (ClimFluids.heat_capacity, state),
                         (ClimFluids.dcons_dq, state), (ClimFluids.modified_chemical_potential, state),
                         (ClimFluids.chemical_potential_derivatives, state_T), (ClimFluids.pressure, state_v),
                         (ClimFluids.chemical_potential, state), (ClimFluids.chemical_potential, state_T))
            @inferred fun(fluid, s)
            @test allocated(fun, fluid, s) == 0
        end
    end

end
