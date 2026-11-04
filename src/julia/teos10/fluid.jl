"""
    fluid = TEOS10Seawater(consvar; F = Float64)

Return an object `fluid` describing seawater with the TEOS-10 thermodynamics, in the float type `F`.
`fluid` can then be used as first argument of thermodynamic functions.
Argument `consvar` defines the conservative variable chosen. Valid values are `:conservative_temperature`,
`:potential_temperature` and `:entropy`.

The thermodynamic potential is the specific enthalpy

    h(p, Θ, q) = Cp0 (Θ - T0) + ∫ v(p', Θ, q) dp'   (from p0 to p)

with v the 75-term polynomial of Roquet et al. (2015), together with the entropy η(Θ, q) and the potential temperature
θ(Θ, q) of the TEOS-10 Gibbs function at the pressure p0 (IAPWS-08, IAPWS-09). All state functions derive from this
potential, so they are thermodynamically consistent; the fluid departs from the full TEOS-10 Gibbs function by the
error of the 75-term fit, valid for 0 ≤ q ≤ 0.042, -2 °C ≤ Θ ≤ 40 °C and p - p0 ≤ 10⁸ Pa.

Units are SI: absolute pressure `p` [Pa], absolute salinity `q` [kg/kg], temperatures `T`, `θ`, `Θ` [K] (Θ in K is the
conservative temperature in °C plus T0 = 273.15 K), entropy [J/kg/K]. The normal pressure is `fluid.p0` = 101325 Pa.
"""
struct TEOS10Seawater{CV, F} <: BinaryFluid{F}
    p0::F       # normal pressure [Pa]
    T0::F       # Celsius zero point [K]
    Cp0::F      # potential enthalpy per unit conservative temperature [J/kg/K]
    q_star::F   # reducing salinity of the Gibbs function [kg/kg]
    function TEOS10Seawater(consvar::Symbol; F = Float64)
        @assert consvar in (:conservative_temperature, :potential_temperature, :entropy)
        return new{consvar, F}(F(101325), F(273.15), F(3991.86795711963), F(35.16504) / 875)
    end
end

const TEOS  = TEOS10Seawater
const TEOSC = TEOS10Seawater{:conservative_temperature}
const TEOSP = TEOS10Seawater{:potential_temperature}
const TEOSS = TEOS10Seawater{:entropy}

@fastmath @muladd @inlineall begin

    ## the conservative variable φ and the temperatures θ, Θ at p0

    # θ(φ, q), Θ(φ, q), ∂Θ/∂φ, ∂Θ/∂q at fixed φ: dh0 = Cp0 dΘ = θ dη + μ0 dq = cp0 dθ + (∂h0/∂q)_θ dq
    TEOS_θΘ_φq(fluid::TEOSC, Θ, q) = (TEOS_θ_Θq(fluid, Θ, q), Θ, one(Θ), zero(Θ))
    function TEOS_θΘ_φq(fluid::TEOSP, θ, q)
        (; T0, Cp0) = fluid
        g, g_T, g_q, g_TT, g_Tq = TEOS_gibbs0(fluid, θ, q)
        return θ, T0 + (g - θ * g_T) / Cp0, -θ * g_TT / Cp0, (g_q - θ * g_Tq) / Cp0
    end
    function TEOS_θΘ_φq(fluid::TEOSS, η, q)
        (; T0, Cp0) = fluid
        θ = TEOS_θ_ηq(fluid, η, q)
        g, g_T, g_q = TEOS_gibbs0(fluid, θ, q)
        return θ, T0 + (g - θ * g_T) / Cp0, θ / Cp0, g_q / Cp0
    end

    # Θ(φ, q), ∂Θ/∂φ, ∂Θ/∂q, without θ when φ = Θ
    TEOS_Θ_φq(fluid::TEOSC, Θ, q) = (Θ, one(Θ), zero(Θ))
    TEOS_Θ_φq(fluid::TEOS, φ, q) = Base.tail(TEOS_θΘ_φq(fluid, φ, q))

    # θ(φ, q)
    TEOS_θ_φq(fluid::TEOS, φ, q)  = TEOS_θΘ_φq(fluid, φ, q)[1]
    TEOS_θ_φq(fluid::TEOSP, θ, q) = θ

    # φ(θ, q)
    TEOS_φ_θq(fluid::TEOSC, θ, q) = TEOS_Θ_θq(fluid, θ, q)
    TEOS_φ_θq(fluid::TEOSP, θ, q) = θ
    TEOS_φ_θq(fluid::TEOSS, θ, q) = TEOS_η_θq(fluid, θ, q)

    # φ(Θ, q)
    TEOS_φ_Θq(fluid::TEOSC, Θ, q) = Θ
    TEOS_φ_Θq(fluid::TEOS, Θ, q)  = TEOS_φ_θq(fluid, TEOS_θ_Θq(fluid, Θ, q), q)

    ## the state at (p, Θ, q), with θ(Θ, q)

    # T(p, θ, Θ, q) = (∂h/∂η)_p,q = (∂h/∂Θ) θ / Cp0, since ∂η/∂Θ = Cp0 / θ
    TEOS_T(fluid, xyz, θ) = TEOS_h_Θ(fluid, xyz) * θ / fluid.Cp0

    # T and its derivatives in (p, Θ, q), with ∂θ/∂Θ = Cp0 / cp0 and ∂θ/∂q = -(∂h0/∂q)_θ / cp0, cp0 = -θ ∂²g/∂T²(θ)
    function TEOS_T_derivatives(fluid, xyz, θ, q)
        (; Cp0) = fluid
        _, g_T, g_q, g_TT, g_Tq = TEOS_gibbs0(fluid, θ, q)
        cp0 = -θ * g_TT
        θ_Θ, θ_q = Cp0 / cp0, -(g_q - θ * g_Tq) / cp0
        h_Θ, h_ΘΘ, h_Θq = TEOS_h_Θ(fluid, xyz), TEOS_h_ΘΘ(fluid, xyz), TEOS_h_Θq(fluid, xyz)
        T = h_Θ * θ / Cp0
        T_p = TEOS_v_Θ(fluid, xyz) * θ / Cp0
        T_Θ = (h_ΘΘ * θ + h_Θ * θ_Θ) / Cp0
        T_q = (h_Θq * θ + h_Θ * θ_q) / Cp0
        return T, T_p, T_Θ, T_q
    end

    ## inversions, by Newton iterations from fixed starting points

    # p(v, Θ, q): v(p, Θ, q) = v from p = p0, for p - p0 up to 10⁸ Pa
    function TEOS_p_vΘq(fluid::AbstractFluid{F}, v, Θ, q) where {F}
        p = fluid.p0 + zero(v * Θ * q)
        for _ in 1:TEOS_newton_steps(F, 2e-9, 1e8, 1e8)
            xyz = TEOS_xyz(fluid, p, Θ, q)
            p = p - (TEOS_v(fluid, xyz) - v) / TEOS_v_p(fluid, xyz)
        end
        return p
    end

    # Θ(p, v, q): v(p, Θ, q) = v from Θ = 25 °C, the root on the warm side of the density maximum; v is convex in Θ,
    # and the first 4 steps bring Θ near the root where ∂v/∂Θ is small (cold water)
    function TEOS_Θ_pvq(fluid::AbstractFluid{F}, p, v, q) where {F}
        Θ = fluid.T0 + 25 + zero(p * v * q)
        for _ in 1:TEOS_newton_steps(F, 0.05, 10, 300) + 4
            xyz = TEOS_xyz(fluid, p, Θ, q)
            Θ = Θ - (TEOS_v(fluid, xyz) - v) / TEOS_v_Θ(fluid, xyz)
        end
        return Θ
    end

    # θ(p, T, q): T(p, Θ(θ, q), q) = T from θ = T, with dT/dθ = (∂T/∂Θ) cp0 / Cp0
    function TEOS_θ_pTq(fluid::AbstractFluid{F}, p, T, q) where {F}
        θ = T + zero(p * q)
        salinity = TEOS_salinity(fluid, q)
        for _ in 1:TEOS_newton_steps(F, 1e-3, 2, 300)
            g, g_T, _, g_TT = TEOS_gibbs0(fluid, θ, salinity)
            xyz = TEOS_xyz(fluid, p, fluid.T0 + (g - θ * g_T) / fluid.Cp0, q)
            h_Θ = TEOS_h_Θ(fluid, xyz)
            T_θ = (TEOS_h_ΘΘ(fluid, xyz) * θ * (-θ * g_TT) / fluid.Cp0 + h_Θ) / fluid.Cp0
            θ = θ - (h_Θ * θ / fluid.Cp0 - T) / T_θ
        end
        return θ
    end

    # (p, θ)(v, T, q): v(p, Θ(θ, q), q) = v and T(p, Θ(θ, q), q) = T from (p0, T), Newton in two dimensions
    function TEOS_pθ_vTq(fluid::AbstractFluid{F}, v, T, q) where {F}
        (; Cp0) = fluid
        p, θ = fluid.p0 + zero(v * T * q), T + zero(v * q)
        salinity = TEOS_salinity(fluid, q)
        for _ in 1:TEOS_newton_steps(F, 2e-9, 1e8, 1e8) + 1
            g, g_T, _, g_TT = TEOS_gibbs0(fluid, θ, salinity)
            Θ_θ = -θ * g_TT / Cp0
            xyz = TEOS_xyz(fluid, p, fluid.T0 + (g - θ * g_T) / Cp0, q)
            h_Θ, v_Θ = TEOS_h_Θ(fluid, xyz), TEOS_v_Θ(fluid, xyz)
            r_v, r_T = TEOS_v(fluid, xyz) - v, h_Θ * θ / Cp0 - T
            J_vp, J_vθ = TEOS_v_p(fluid, xyz), v_Θ * Θ_θ
            J_Tp, J_Tθ = v_Θ * θ / Cp0, (TEOS_h_ΘΘ(fluid, xyz) * Θ_θ * θ + h_Θ) / Cp0
            det = J_vp * J_Tθ - J_vθ * J_Tp
            p = p - (J_Tθ * r_v - J_vθ * r_T) / det
            θ = θ - (J_vp * r_T - J_Tp * r_v) / det
        end
        return p, θ
    end

    ## state functions at (p, φ, q)

    conservative_variable(      fluid::TEOS, (p, φ, q)::PConsQ) = φ
    pressure(                   fluid::TEOS, (p, φ, q)::PConsQ) = p
    potential_temperature(      fluid::TEOS, (p, φ, q)::PConsQ) = TEOS_θ_φq(fluid, φ, q)
    specific_entropy(           fluid::TEOS, (p, φ, q)::PConsQ) = TEOS_η_θq(fluid, TEOS_θ_φq(fluid, φ, q), q)
    specific_entropy(           fluid::TEOSS, (p, η, q)::PConsQ) = η

    function specific_enthalpy(fluid::TEOS, (p, φ, q)::PConsQ)
        Θ = TEOS_Θ_φq(fluid, φ, q)[1]
        return TEOS_h(fluid, TEOS_xyz(fluid, p, Θ, q), Θ)
    end
    specific_volume(fluid::TEOS, (p, φ, q)::PConsQ) = TEOS_v(fluid, TEOS_xyz(fluid, p, TEOS_Θ_φq(fluid, φ, q)[1], q))
    potential_volume(fluid::TEOS, (p, φ, q)::PConsQ) = specific_volume(fluid, (; p = fluid.p0, consvar = φ, q))
    potential_enthalpy(fluid::TEOS, (p, φ, q)::PConsQ) = fluid.Cp0 * (TEOS_Θ_φq(fluid, φ, q)[1] - fluid.T0)
    function specific_internal_energy(fluid::TEOS, (p, φ, q)::PConsQ)
        Θ = TEOS_Θ_φq(fluid, φ, q)[1]
        xyz = TEOS_xyz(fluid, p, Θ, q)
        return TEOS_h(fluid, xyz, Θ) - p * TEOS_v(fluid, xyz)
    end
    # c² = -v² / (∂v/∂p)_η,q
    function sound_speed2(fluid::TEOS, (p, φ, q)::PConsQ)
        xyz = TEOS_xyz(fluid, p, TEOS_Θ_φq(fluid, φ, q)[1], q)
        return -TEOS_v(fluid, xyz)^2 / TEOS_v_p(fluid, xyz)
    end

    # the conjugate variable ∂h/∂φ = (∂h/∂Θ) ∂Θ/∂φ, which is T for φ = η
    function conjugate_variable(fluid::TEOS, (p, φ, q)::PConsQ)
        Θ, Θ_φ = TEOS_Θ_φq(fluid, φ, q)
        return TEOS_h_Θ(fluid, TEOS_xyz(fluid, p, Θ, q)) * Θ_φ
    end
    # (∂h/∂q)_p,φ = (∂h/∂q)_Θ + (∂h/∂Θ) (∂Θ/∂q)_φ
    function modified_chemical_potential(fluid::TEOS, (p, φ, q)::PConsQ)
        Θ, _, Θ_q = TEOS_Θ_φq(fluid, φ, q)
        xyz = TEOS_xyz(fluid, p, Θ, q)
        return TEOS_h_q(fluid, xyz) + TEOS_h_Θ(fluid, xyz) * Θ_q
    end
    # h, v, ∂h/∂φ
    function exner_functions(fluid::TEOS, (p, φ, q)::PConsQ)
        Θ, Θ_φ = TEOS_Θ_φq(fluid, φ, q)
        xyz = TEOS_xyz(fluid, p, Θ, q)
        return TEOS_h(fluid, xyz, Θ), TEOS_v(fluid, xyz), TEOS_h_Θ(fluid, xyz) * Θ_φ
    end
    # v, ∂v/∂p, ∂v/∂φ, ∂v/∂q (at fixed φ)
    function volume_functions(fluid::TEOS, (p, φ, q)::PConsQ)
        Θ, Θ_φ, Θ_q = TEOS_Θ_φq(fluid, φ, q)
        xyz = TEOS_xyz(fluid, p, Θ, q)
        v_Θ = TEOS_v_Θ(fluid, xyz)
        return TEOS_v(fluid, xyz), TEOS_v_p(fluid, xyz), v_Θ * Θ_φ, TEOS_v_q(fluid, xyz) + v_Θ * Θ_q
    end

    function temperature(fluid::TEOS, (p, φ, q)::PConsQ)
        θ, Θ = TEOS_θΘ_φq(fluid, φ, q)
        return TEOS_T(fluid, TEOS_xyz(fluid, p, Θ, q), θ)
    end
    # cp = T (∂η/∂T)_p,q = (∂h/∂Θ) / (∂T/∂Θ)
    function heat_capacity(fluid::TEOS, (p, φ, q)::PConsQ)
        θ, Θ = TEOS_θΘ_φq(fluid, φ, q)
        xyz = TEOS_xyz(fluid, p, Θ, q)
        _, _, T_Θ = TEOS_T_derivatives(fluid, xyz, θ, q)
        return TEOS_h_Θ(fluid, xyz) / T_Θ
    end
    # (∂φ/∂q)_p,T = ((∂Θ/∂q)_p,T - (∂Θ/∂q)_φ) / (∂Θ/∂φ), with (∂Θ/∂q)_p,T = -(∂T/∂q)_Θ / (∂T/∂Θ)
    function dcons_dq(fluid::TEOS, (p, φ, q)::PConsQ)
        θ, Θ, Θ_φ, Θ_q = TEOS_θΘ_φq(fluid, φ, q)
        xyz = TEOS_xyz(fluid, p, Θ, q)
        _, _, T_Θ, T_q = TEOS_T_derivatives(fluid, xyz, θ, q)
        return (-T_q / T_Θ - Θ_q) / Θ_φ
    end

    # μ = (∂g/∂q)_p,T = (∂h/∂q)_p,η = (∂h/∂q)_Θ + (∂h/∂Θ) μ0 / Cp0, with μ0 = ∂g/∂q(θ) at p0 [J/kg]
    TEOS_μ(fluid, xyz, g_q) = TEOS_h_q(fluid, xyz) + TEOS_h_Θ(fluid, xyz) * g_q / fluid.Cp0
    TEOS_μ(fluid, p, θ, Θ, q) = TEOS_μ(fluid, TEOS_xyz(fluid, p, Θ, q), TEOS_gibbs0(fluid, θ, q)[3])

    function chemical_potential(fluid::TEOS, (p, φ, q)::PConsQ)
        θ, Θ = TEOS_θΘ_φq(fluid, φ, q)
        return TEOS_μ(fluid, p, θ, Θ, q)
    end
    function chemical_potential(fluid::TEOS, (p, T, q)::PTQ)
        θ = TEOS_θ_pTq(fluid, p, T, q)
        return TEOS_μ(fluid, p, θ, TEOS_Θ_θq(fluid, θ, q), q)
    end

    # μ and its derivatives at fixed (p, T)
    function chemical_potential_derivatives(fluid::TEOS, (p, φ, q)::PConsQ)
        θ, Θ = TEOS_θΘ_φq(fluid, φ, q)
        return TEOS_μ_derivatives(fluid, p, θ, Θ, q)
    end
    function chemical_potential_derivatives(fluid::TEOS, (p, T, q)::PTQ)
        θ = TEOS_θ_pTq(fluid, p, T, q)
        return TEOS_μ_derivatives(fluid, p, θ, TEOS_Θ_θq(fluid, θ, q), q)
    end
    function TEOS_μ_derivatives(fluid, p, θ, Θ, q)
        (; Cp0) = fluid
        xyz = TEOS_xyz(fluid, p, Θ, q)
        _, g_T, g_q, g_TT, g_Tq, g_qq = TEOS_gibbs0(fluid, θ, q)
        cp0 = -θ * g_TT
        θ_Θ, θ_q = Cp0 / cp0, -(g_q - θ * g_Tq) / cp0
        h_Θ, h_ΘΘ, h_Θq = TEOS_h_Θ(fluid, xyz), TEOS_h_ΘΘ(fluid, xyz), TEOS_h_Θq(fluid, xyz)
        v_Θ = TEOS_v_Θ(fluid, xyz)
        # μ and its derivatives in (p, Θ, q)
        μ   = TEOS_μ(fluid, xyz, g_q)
        μ_p = TEOS_v_q(fluid, xyz) + v_Θ * g_q / Cp0
        μ_Θ = h_Θq + (h_ΘΘ * g_q + h_Θ * g_Tq * θ_Θ) / Cp0
        μ_q = TEOS_h_qq(fluid, xyz) + (h_Θq * g_q + h_Θ * (g_qq + g_Tq * θ_q)) / Cp0
        # T(p, Θ, q) and its derivatives, to change to (p, T, q)
        T_p = v_Θ * θ / Cp0
        T_Θ = (h_ΘΘ * θ + h_Θ * θ_Θ) / Cp0
        T_q = (h_Θq * θ + h_Θ * θ_q) / Cp0
        return μ, μ_p - μ_Θ * T_p / T_Θ, μ_Θ / T_Θ, μ_q - μ_Θ * T_q / T_Θ
    end

    ## other state variables, converted to (p, φ, q)

    canonical_state(fluid::TEOS, (p, φ, q)::PConsQ)  = (; p, consvar = φ, q)
    canonical_state(fluid::TEOS, (p, T, q)::PTQ)     = (; p, consvar = TEOS_φ_θq(fluid, TEOS_θ_pTq(fluid, p, T, q), q), q)
    canonical_state(fluid::TEOS, (p, θ, q)::PThQ)    = (; p, consvar = TEOS_φ_θq(fluid, θ, q), q)
    canonical_state(fluid::TEOS, (p, s, q)::PSQ)     = (; p, consvar = TEOS_φ_θq(fluid, TEOS_θ_ηq(fluid, s, q), q), q)
    canonical_state(fluid::TEOS, (p, v, q)::PVQ)     = (; p, consvar = TEOS_φ_Θq(fluid, TEOS_Θ_pvq(fluid, p, v, q), q), q)
    canonical_state(fluid::TEOS, (v, φ, q)::VConsQ)  = (; p = TEOS_p_vΘq(fluid, v, TEOS_Θ_φq(fluid, φ, q)[1], q), consvar = φ, q)
    function canonical_state(fluid::TEOS, (v, s, q)::VSQ)
        θ = TEOS_θ_ηq(fluid, s, q)
        return (; p = TEOS_p_vΘq(fluid, v, TEOS_Θ_θq(fluid, θ, q), q), consvar = TEOS_φ_θq(fluid, θ, q), q)
    end
    function canonical_state(fluid::TEOS, (v, T, q)::VTQ)
        p, θ = TEOS_pθ_vTq(fluid, v, T, q)
        return (; p, consvar = TEOS_φ_θq(fluid, θ, q), q)
    end
    pressure(fluid::TEOS, (v, φ, q)::VConsQ) = TEOS_p_vΘq(fluid, v, TEOS_Θ_φq(fluid, φ, q)[1], q)

end
