# The Gibbs function g(q, T, p0) of seawater at the normal pressure and its derivatives, with q the absolute salinity
# [kg/kg] and T the temperature [K]. With the polynomials cᵢ(τ) of `coefficients.jl` and a(τ) = g₁₀₀ + g₁₁₀ τ,
#   g     = gW(τ) + a ξ² ln ξ + Σᵢ cᵢ ξⁱ
#   ∂g/∂q = (a (2 ln ξ + 1) + Σᵢ i cᵢ ξⁱ⁻²) / 2q*,          ∂²g/∂q² = (2a + Σᵢ i (i - 2) cᵢ ξⁱ⁻²) / 4q*²ξ²
# with ∂/∂T = (1/T*) ∂/∂τ. The logarithm is evaluated at ξ² + eps, so that g and ∂g/∂T are finite at q = 0.

@fastmath @muladd @inlineall begin

    # evalpoly with the coefficients converted to the float type F of the fluid
    TEOS_poly(::Type{F}, x, c) where {F} = evalpoly(x, map(F, c))

    # p(τ), p'(τ), p''(τ) from the coefficients of p, p', p''
    TEOS_poly_derivatives(::Type{F}, τ, (c, c_τ, c_ττ)) where {F} = (TEOS_poly(F, τ, c), TEOS_poly(F, τ, c_τ), TEOS_poly(F, τ, c_ττ))

    # (ξ², ξ, ln ξ), which depend on q only: computed once in iterations over T
    function TEOS_salinity(fluid::AbstractFluid{F}, q) where {F}
        ξ2 = q / fluid.q_star
        return ξ2, sqrt(ξ2), log(ξ2 + eps(F)) / 2
    end

    # g, ∂g/∂T, ∂g/∂q, ∂²g/∂T², ∂²g/∂T∂q, ∂²g/∂q² at (q, T, p0)
    TEOS_gibbs0(fluid, T, q) = TEOS_gibbs0(fluid, T, TEOS_salinity(fluid, q))
    function TEOS_gibbs0(fluid::AbstractFluid{F}, T, (ξ2, ξ, L)::Tuple) where {F}
        (; T0, q_star) = fluid
        T_star = F(40)
        τ = (T - T0) / T_star

        # pure water
        gW, gW_τ, gW_ττ = TEOS_poly_derivatives(F, τ, TEOS_GW_ALL)

        # saline part: G = a ξ² ln ξ + Σᵢ cᵢ ξⁱ, written as ξ² (a ln ξ + P(ξ)), P = Σᵢ cᵢ ξⁱ⁻²,
        # and the sums Σᵢ i cᵢ ξⁱ⁻² (for ∂/∂q) and Σᵢ i (i - 2) cᵢ ξⁱ⁻² (for ∂²/∂q²)
        a, a_τ = F(TEOS_GS_LOG[1]) + F(TEOS_GS_LOG[2]) * τ, F(TEOS_GS_LOG[2])
        (c2, c2_τ, c2_ττ), (c3, c3_τ, c3_ττ), (c4, c4_τ, c4_ττ), (c5, c5_τ, c5_ττ), (c6, c6_τ, c6_ττ), (c7, c7_τ, c7_ττ) =
            map(c -> TEOS_poly_derivatives(F, τ, c), TEOS_GS_ALL)
        P    = evalpoly(ξ, (c2, c3, c4, c5, c6, c7))
        P_τ  = evalpoly(ξ, (c2_τ, c3_τ, c4_τ, c5_τ, c6_τ, c7_τ))
        P_ττ = evalpoly(ξ, (c2_ττ, c3_ττ, c4_ττ, c5_ττ, c6_ττ, c7_ττ))
        Q    = evalpoly(ξ, (2c2, 3c3, 4c4, 5c5, 6c6, 7c7))
        Q_τ  = evalpoly(ξ, (2c2_τ, 3c3_τ, 4c4_τ, 5c5_τ, 6c6_τ, 7c7_τ))
        R    = evalpoly(ξ, (3c3, 8c4, 15c5, 24c6, 35c7))         # Σᵢ i (i - 2) cᵢ ξⁱ⁻³

        g     = gW + ξ2 * (a * L + P)
        g_T   = (gW_τ + ξ2 * (a_τ * L + P_τ)) / T_star
        g_TT  = (gW_ττ + ξ2 * P_ττ) / T_star^2
        g_q   = (a * (2L + 1) + Q) / (2q_star)
        g_Tq  = (a_τ * (2L + 1) + Q_τ) / (2q_star * T_star)
        g_qq  = (2a + ξ * R) / (4q_star^2 * ξ2)
        return g, g_T, g_q, g_TT, g_Tq, g_qq
    end

end
