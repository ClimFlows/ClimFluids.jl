# Time per point of the state functions of `TEOS10Seawater`, for its three conservative variables, and of
# `NonlinearBinaryFluid` for comparison. Needs ClimFluids and BenchmarkTools in the active environment:
#   julia -O3 benchmark/teos10.jl
using ClimFluids, BenchmarkTools, Printf

# the thermodynamics of a cell of a compressible model, from (v, φ, q)
function cell_thermodynamics(fluid, v, φ, q)
    p = fluid(:v, :consvar, :q).pressure(v, φ, q)
    fluid_pφq = fluid(:p, :consvar, :q)
    T = fluid_pφq.temperature(p, φ, q)
    _, μ_p, _, μ_q = fluid(:p, :T, :q).chemical_potential_derivatives(p, T, q)
    return p + T + fluid_pφq.conjugate_variable(p, φ, q) + fluid_pφq.modified_chemical_potential(p, φ, q) +
           μ_p / μ_q + fluid_pφq.heat_capacity(p, φ, q) + fluid_pφq.dcons_dq(p, φ, q) + fluid_pφq.specific_enthalpy(p, φ, q)
end

specific_volume(fluid, p, φ, q) = fluid(:p, :consvar, :q).specific_volume(p, φ, q)
volume_functions(fluid, p, φ, q) = sum(fluid(:p, :consvar, :q).volume_functions(p, φ, q))
exner_functions(fluid, p, φ, q) = sum(fluid(:p, :consvar, :q).exner_functions(p, φ, q))
temperature(fluid, p, φ, q) = fluid(:p, :consvar, :q).temperature(p, φ, q)
pressure(fluid, v, φ, q) = fluid(:v, :consvar, :q).pressure(v, φ, q)

function apply!(out, f::G, fluid, a, b, c) where {G}
    @inbounds @simd for i in eachindex(out, a, b, c)
        out[i] = f(fluid, a[i], b[i], c[i])
    end
    return out
end

bytes(args...) = (apply!(args...); @allocated apply!(args...))

function benchmark(fluid; N = 10_000)
    p = 101325 .+ 1e6 * rand(N)                    # down to 100 m
    q = 0.03 .+ 0.008 * rand(N)
    φ = [fluid(:p, :T, :q).conservative_variable(p[i], 275 + 25rand(), q[i]) for i in 1:N]
    v = [fluid(:p, :consvar, :q).specific_volume(p[i], φ[i], q[i]) for i in 1:N]
    out = zeros(N)
    println(typeof(fluid))
    for (f, a) in ((specific_volume, p), (volume_functions, p), (exner_functions, p), (temperature, p),
                   (pressure, v), (cell_thermodynamics, v))
        t = @belapsed apply!($out, $f, $fluid, $a, $φ, $q)
        @printf "  %-20s %7.1f ns/point, %d bytes\n" f t / N * 1e9 bytes(out, f, fluid, a, φ, q)
    end
end

for consvar in (:conservative_temperature, :potential_temperature, :entropy)
    benchmark(TEOS10Seawater(consvar))
end
benchmark(NonlinearBinaryFluid(:potential_temperature, (; p0 = 1e5, T0 = 283.0, q0 = 0.035, Cp0 = 4000.0, v0 = 1e-3,
    α_T = 2e-4, α_q = 0.8, α_p = 4.4e-10, α_TT = 1e-5, γ = 1e-8, M_s = 3.14038218e-2, R = 8.31446261815324)))
