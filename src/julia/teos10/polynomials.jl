# The dynamic enthalpy H(x, y, z) = ∫₀ᶻ v dz' = Σ Vᵢⱼₖ xⁱ yʲ zᵏ⁺¹ / (k + 1), so that the specific enthalpy is
#   h(p, Θ, q) = Cp0 (Θ - T0) + 10⁸ Pa × H,   and v = ∂h/∂p = ∂H/∂z,
# to the rounding of the coefficients: all derivatives come from the same coefficients Vᵢⱼₖ / (k + 1) of H.
# `TEOS_dH(F, x, y, z, Val((a, b, c)))` evaluates ∂ᵃ⁺ᵇ⁺ᶜH/∂xᵃ∂yᵇ∂zᶜ by Horner's scheme (z outermost, y innermost), its
# coefficients computed once, when the method is compiled, and converted to the float type F.

# n (n - 1) ⋯ (n - m + 1)
TEOS_falling(n, m) = prod(n-m+1:n; init = 1)

# Horner expression of Σ cᵢⱼₖ xⁱ yʲ zᵏ, with the coefficients in a dictionary (i, j, k) => cᵢⱼₖ
function TEOS_horner(F, coefs)
    horner(var, terms) = foldr((t, acc) -> :(muladd($acc, $var, $t)), terms[1:end-1]; init = terms[end])
    degree(d) = maximum(key -> key[d], keys(coefs); init = 0)
    poly_y(i, k) = horner(:y, [F(get(coefs, (i, j, k), 0.0)) for j in 0:degree(2)])
    poly_x(k) = horner(:x, [poly_y(i, k) for i in 0:degree(1)])
    return horner(:z, [poly_x(k) for k in 0:degree(3)])
end

@generated function TEOS_dH(::Type{F}, x, y, z, ::Val{abc}) where {F, abc}
    a, b, c = abc
    coefs = Dict{NTuple{3, Int}, Float64}()
    for (i, j, k, V) in TEOS_V
        if i >= a && j >= b && k + 1 >= c
            coefs[(i - a, j - b, k + 1 - c)] = V / (k + 1) * TEOS_falling(i, a) * TEOS_falling(j, b) * TEOS_falling(k + 1, c)
        end
    end
    isempty(coefs) && return :(zero(x * y * z))
    return Expr(:block, Expr(:meta, :inline), TEOS_horner(F, coefs))
end

@fastmath @muladd @inlineall begin

    # (x, y, z, ∂x/∂q) of the polynomial at (p, Θ, q); the methods below take this tuple, `xyz`
    function TEOS_xyz(fluid::AbstractFluid{F}, p, Θ, q) where {F}
        (; p0, T0, q_star) = fluid
        x = sqrt((q + F(0.024)) / q_star)
        return x, (Θ - T0) / 40, (p - p0) / 100_000_000, 1 / (2q_star * x)
    end

    # derivatives of h(p, Θ, q), from those of H
    TEOS_H(fluid::AbstractFluid{F}, (x, y, z, _), a, b, c) where {F} = TEOS_dH(F, x, y, z, Val((a, b, c)))

    # h, v
    TEOS_h(fluid, xyz, Θ)   = fluid.Cp0 * (Θ - fluid.T0) + 100_000_000 * TEOS_H(fluid, xyz, 0, 0, 0)
    TEOS_v(fluid, xyz)      = TEOS_H(fluid, xyz, 0, 0, 1)
    # ∂h/∂Θ, ∂h/∂q
    TEOS_h_Θ(fluid, xyz)    = fluid.Cp0 + 2_500_000 * TEOS_H(fluid, xyz, 0, 1, 0)
    TEOS_h_q(fluid, xyz)    = 100_000_000 * xyz[4] * TEOS_H(fluid, xyz, 1, 0, 0)
    # ∂v/∂p, ∂v/∂Θ, ∂v/∂q
    TEOS_v_p(fluid, xyz)    = TEOS_H(fluid, xyz, 0, 0, 2) / 100_000_000
    TEOS_v_Θ(fluid, xyz)    = TEOS_H(fluid, xyz, 0, 1, 1) / 40
    TEOS_v_q(fluid, xyz)    = xyz[4] * TEOS_H(fluid, xyz, 1, 0, 1)
    # ∂²h/∂Θ², ∂²h/∂Θ∂q, ∂²h/∂q², with ∂²x/∂q² = -(∂x/∂q)² / x
    TEOS_h_ΘΘ(fluid, xyz)   = 62_500 * TEOS_H(fluid, xyz, 0, 2, 0)
    TEOS_h_Θq(fluid, xyz)   = 2_500_000 * xyz[4] * TEOS_H(fluid, xyz, 1, 1, 0)
    TEOS_h_qq(fluid, xyz)   = 100_000_000 * xyz[4]^2 * (TEOS_H(fluid, xyz, 2, 0, 0) - TEOS_H(fluid, xyz, 1, 0, 0) / xyz[1])

end
