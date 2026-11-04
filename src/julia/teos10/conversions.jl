# Conversions between potential temperature θ [K], conservative temperature Θ [K] and entropy η [J/kg/K] at fixed
# absolute salinity q [kg/kg], all at the normal pressure p0, from the Gibbs function g(q, T, p0):
#   η = -∂g/∂T(q, θ, p0),    Cp0 (Θ - T0) = h0 = g - θ ∂g/∂T  (the potential enthalpy),
# with ∂h0/∂θ = θ ∂η/∂θ = -θ ∂²g/∂T² = cp(q, θ, p0). The inverses θ(Θ, q) and θ(η, q) are found by Newton iterations
# from accurate initial guesses, with a fixed number of steps so that they are free of branches.

@fastmath @muladd @inlineall begin

    # The number of Newton steps that bring an error e₀ down to the round-off of the float type F for a variable of
    # magnitude X, if each step maps the error e to C e²: after n steps, C e = (C e₀)^(2ⁿ).
    TEOS_newton_steps(::Type{F}, C, e₀, X) where {F} = ceil(Int, log2(log(C * eps(F) * X) / log(C * e₀)))

    # h0(θ, q), the potential enthalpy [J/kg]
    function TEOS_h0_θq(fluid, θ, q)
        g, g_T = TEOS_gibbs0(fluid, θ, q)
        return g - θ * g_T
    end
    # Θ(θ, q)
    TEOS_Θ_θq(fluid, θ, q) = fluid.T0 + TEOS_h0_θq(fluid, θ, q) / fluid.Cp0
    # η(θ, q)
    TEOS_η_θq(fluid, θ, q) = -TEOS_gibbs0(fluid, θ, q)[2]

    # θ(Θ, q): Newton on h0(θ, q) = Cp0 (Θ - T0), from a rational function of Θ and q (error below 0.01 K)
    function TEOS_θ_Θq(fluid::AbstractFluid{F}, Θ, q) where {F}
        (; T0, Cp0) = fluid
        s = 1000q / F(1.0047154285714286)                                 # S_A / (35.16504/35 g/kg)
        t = Θ - T0                                                        # Θ [°C]
        num = F(-1.446013646344788e-2) + s * (F(-3.305308995852924e-3) + F(1.062415929128982e-4) * s) +
              t * (F(9.477566673794488e-1) + F(2.166591947736613e-3) * s + F(3.828842955039902e-3) * t)
        den = 1 + F(6.506097115635800e-4) * s + t * (F(3.830289486850898e-3) + F(1.247811760368034e-6) * t)
        θ = T0 + num / den
        h0 = Cp0 * t
        salinity = TEOS_salinity(fluid, q)
        for _ in 1:TEOS_newton_steps(F, 1e-3, 0.02, 300)
            g, g_T, _, g_TT = TEOS_gibbs0(fluid, θ, salinity)
            θ = θ + (g - θ * g_T - h0) / (θ * g_TT)
        end
        return θ
    end

    # θ(η, q): Newton on η(θ, q) = η, from a least-squares fit of θ(η, q) over the range of the 75-term polynomial,
    #   θ ≈ T0 + Σᵢⱼ cᵢⱼ eⁱ ξʲ + (c₀ + c₁ e) ξ² ln ξ,   e = η / (600 J/kg/K),
    # with an error below 0.01 K (0.03 K up to 45 °C and 50 g/kg)
    function TEOS_θ_ηq(fluid::AbstractFluid{F}, η, q) where {F}
        salinity = TEOS_salinity(fluid, q)
        ξ2, ξ, L = salinity
        e = η / 600
        c0 = evalpoly(ξ, (F(0.00959935509508146), zero(F), F(0.1703852080031997), F(-0.5256059390678489), F(0.46002731008497677)))
        c1 = evalpoly(ξ, (F(38.864791433115975), zero(F), F(2.610000798148543), F(-0.10405556419435828)))
        c2 = evalpoly(ξ, (F(3.2402482178289307), zero(F), F(-0.03518157606259137)))
        θ = fluid.T0 + evalpoly(e, (c0, c1, c2, F(-0.03843044750232895))) +
            ξ2 * L * (F(1.3420739401078954) + F(0.1796601908462941) * e)
        for _ in 1:TEOS_newton_steps(F, 2e-3, 0.01, 300)
            _, g_T, _, g_TT = TEOS_gibbs0(fluid, θ, salinity)
            θ = θ - (g_T + η) / g_TT
        end
        return θ
    end

end
