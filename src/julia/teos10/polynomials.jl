# The dynamic enthalpy H(x, y, z) = ∫₀ᶻ v dz' = Σ Hᵢⱼₖ xⁱ yʲ zᵏ⁺¹, with Hᵢⱼₖ = Vᵢⱼₖ / (k + 1), gives
#   h(p, Θ, q) = Cp0 (Θ - T0) + 10⁸ Pa × H   and   v = ∂h/∂p = ∂H/∂z.
# `TEOS_dH(V, H, x, y, z, Val((a, b, d)))` evaluates ∂ᵃ⁺ᵇ⁺ᵈH/∂xᵃ∂yᵇ∂zᵈ by Horner's scheme (z outermost, y innermost),
# from the coefficients Vᵢⱼₖ and Hᵢⱼₖ held by the fluid.

# n (n - 1) ⋯ (n - m + 1)
TEOS_falling(n, m) = prod(n-m+1:n; init = 1)

# Horner expression of Σ cᵢⱼₖ xⁱ yʲ zᵏ, from a dictionary (i, j, k) => expression of cᵢⱼₖ
function TEOS_horner(zero_coef, coefs)
    horner(var, terms) = foldr((t, acc) -> :(muladd($acc, $var, $t)), terms[1:end-1]; init = terms[end])
    degree(d) = maximum(key -> key[d], keys(coefs); init = 0)
    poly_y(i, k) = horner(:y, [get(coefs, (i, j, k), zero_coef) for j in 0:degree(2)])
    poly_x(k) = horner(:x, [poly_y(i, k) for i in 0:degree(1)])
    return horner(:z, [poly_x(k) for k in 0:degree(3)])
end

@generated function TEOS_dH(V::NTuple{N, F}, H::NTuple{N, F}, x, y, z, ::Val{abc}) where {N, F, abc}
    a, b, d = abc
    coefs = Dict{NTuple{3, Int}, Any}()
    for (n, (i, j, k, _)) in enumerate(TEOS_V)
        if i >= a && j >= b && k + 1 >= d
            # coefficient of xⁱ⁻ᵃ yʲ⁻ᵇ zᵏ⁺¹⁻ᵈ: Hᵢⱼₖ, or Vᵢⱼₖ = (k + 1) Hᵢⱼₖ if d ≥ 1, times integer factors
            factor = TEOS_falling(i, a) * TEOS_falling(j, b) * (d == 0 ? 1 : TEOS_falling(k, d - 1))
            coef = d == 0 ? :(H[$n]) : :(V[$n])
            coefs[(i - a, j - b, k + 1 - d)] = factor == 1 ? coef : :($factor * $coef)
        end
    end
    isempty(coefs) && return :(zero(x * y * z))
    return Expr(:block, Expr(:meta, :inline), TEOS_horner(:(zero($F)), coefs))
end

@fastmath @muladd @inlineall begin

    # (x, y, z, ∂x/∂q) of the polynomial at (p, Θ, q); the methods below take this tuple, `xyz`
    function TEOS_xyz(fluid::AbstractFluid{F}, p, Θ, q) where {F}
        (; p0, T0, q_star) = fluid
        x = sqrt((q + F(3) / 125) / q_star)               # x = √((q + 0.024) / q*)
        return x, (Θ - T0) / 40, (p - p0) / 100_000_000, 1 / (2q_star * x)
    end

    # derivatives of h(p, Θ, q), from those of H
    TEOS_H(fluid, (x, y, z, _), a, b, c) = TEOS_dH(fluid.V, fluid.H, x, y, z, Val((a, b, c)))

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
