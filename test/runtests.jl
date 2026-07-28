using Test
using LinearAlgebra
using SpheroidalWaveFunctions

const SWF = SpheroidalWaveFunctions

function tuple_isapprox(a, b; rtol=1e-10, atol=1e-12)
    return isapprox(a[1], b[1]; rtol=rtol, atol=atol) && isapprox(a[2], b[2]; rtol=rtol, atol=atol)
end

function relative_residual(terms...)
    T = typeof(float(real(first(terms))))
    return abs(sum(terms)) / (sum(abs, terms) + one(T))
end

function second_derivative_from_first(f, x)
    T = typeof(float(real(x)))
    h = cbrt(eps(T)) * max(one(T), abs(x))

    central(h) = begin
        _, df_plus = f(x + h)
        _, df_minus = f(x - h)
        (df_plus - df_minus) / (2h)
    end

    d2_h = central(h)
    d2_h2 = central(h / 2)
    d2 = (4*d2_h2 - d2_h) / 3

    return d2, abs(d2_h2 - d2_h)
end

function derivative_relative_error(f, x)
    T = typeof(float(real(x)))
    h = cbrt(eps(T)) * max(one(T), abs(x))

    central(h) = begin
        y_plus, _ = f(x + h)
        y_minus, _ = f(x - h)
        (y_plus - y_minus) / (2h)
    end

    dy_h = central(h)
    dy_h2 = central(h / 2)
    dy_fd = (4*dy_h2 - dy_h) / 3
    _, dy = f(x)

    return abs(dy - dy_fd) / (abs(dy_fd) + abs(dy) + one(T))
end

function angular_ode_relative_residual(f, m, lam, gamma2, x)
    y, dy = f(x)
    d2y, d2y_error = second_derivative_from_first(f, x)
    potential = lam - gamma2 * x^2 - m^2 / (1 - x^2)
    t1 = (1 - x^2) * d2y
    t2 = -2*x * dy
    t3 = potential * y
    residual = relative_residual(t1, t2, t3)
    d2y_reliability = d2y_error / (abs(d2y) + one(typeof(float(real(x)))))
    return residual, d2y_reliability
end

function radial_ode_relative_residual(f, m, lam, c, xi)
    y, dy = f(xi)
    d2y, d2y_error = second_derivative_from_first(f, xi)
    potential = c^2 * xi^2 - lam - m^2 / (xi^2 - 1)
    t1 = (xi^2 - 1) * d2y
    t2 = 2*xi * dy
    t3 = potential * y
    residual = relative_residual(t1, t2, t3)
    d2y_reliability = d2y_error / (abs(d2y) + one(typeof(float(real(xi)))))
    return residual, d2y_reliability
end

function oblate_radial_ode_relative_residual(f, m, lam, c, xi)
    y, dy = f(xi)
    d2y, d2y_error = second_derivative_from_first(f, xi)
    potential = c^2 * xi^2 - lam + m^2 / (1 + xi^2)
    t1 = (1 + xi^2) * d2y
    t2 = 2*xi * dy
    t3 = potential * y
    residual = relative_residual(t1, t2, t3)
    d2y_reliability = d2y_error / (abs(d2y) + one(typeof(float(real(xi)))))
    return residual, d2y_reliability
end

function coefficient_recurrence_relative_residual(m, n, c, lam, dr)
    gamma2 = real(c^2)
    start = iseven(n - m) ? 0 : 1
    max_residual = 0.0

    for i in 2:(length(dr) - 1)
        r = start + 2 * (i - 1)
        lower = SWF.γᵣ(m, r, gamma2) * dr[i - 1]
        diagonal = (SWF.βᵣ(m, r, gamma2) - lam) * dr[i]
        upper = SWF.αᵣ(m, r, gamma2) * dr[i + 1]
        residual = lower + diagonal + upper
        scale = abs(lower) + abs(diagonal) + abs(upper) + one(typeof(float(real(residual))))
        max_residual = max(max_residual, abs(residual) / scale)
    end

    return max_residual
end

function scaled_prolate_radial_wronskian(m, n, c, xi)
    y1, dy1 = prolate_radial1(m, n, c, xi)
    y2, dy2 = prolate_radial2(m, n, c, xi)
    return (xi^2 - 1) * (y1 * dy2 - dy1 * y2)
end

function gausslegendre_nodes_weights(n)
    beta = [k / sqrt(4k^2 - 1) for k in 1:(n - 1)]
    values = eigen(SymTridiagonal(zeros(n), beta))
    nodes = values.values
    weights = 2 .* abs2.(values.vectors[1, :])
    return nodes, weights
end

function normalized_inner_product(f, g; nquad=160)
    nodes, weights = gausslegendre_nodes_weights(nquad)
    inner = sum(w * first(f(x)) * first(g(x)) for (x, w) in zip(nodes, weights))
    normf2 = sum(w * abs2(first(f(x))) for (x, w) in zip(nodes, weights))
    normg2 = sum(w * abs2(first(g(x))) for (x, w) in zip(nodes, weights))
    return abs(inner) / sqrt(normf2 * normg2)
end

function angular_c0_constant_mode_is_constant()
    vals = [first(prolate_angular_leg(0, 0, 0.0, x)) for x in (-0.8, -0.2, 0.3, 0.9)]
    return all(v -> isapprox(v, vals[1]; rtol=1e-13, atol=1e-13), vals)
end

function angular_c0_linear_mode_is_linear()
    ratios = [first(prolate_angular_leg(0, 1, 0.0, x)) / x for x in (-0.8, -0.3, 0.2, 0.7)]
    return all(r -> isapprox(r, ratios[1]; rtol=1e-12, atol=1e-12), ratios)
end

@testset "Characteristic values: intrinsic checks" begin
    for (m, n) in ((0, 0), (0, 3), (1, 4), (2, 5))
        @test prolate_cv(m, n, 0.0) == n * (n + 1)
        @test oblate_cv(m, n, 0.0) == n * (n + 1)
    end

    for (m, nmax, c) in ((0, 6, 0.4), (1, 6, 1.2), (2, 7, 2.0))
        seq = prolate_cv_seq(m, nmax, c)
        @test length(seq) == nmax - m + 1
        @test all(isfinite, seq)
        @test issorted(seq)
        @test all(diff(seq) .> 0)
        for (i, n) in enumerate(m:nmax)
            @test isapprox(seq[i], prolate_cv(m, n, c); rtol=1e-9, atol=1e-10)
        end
    end

    for (m, nmax, c) in ((0, 6, 0.4), (1, 6, 1.2), (2, 7, 2.0))
        seq = oblate_cv_seq(m, nmax, c)
        @test length(seq) == nmax - m + 1
        @test all(isfinite, seq)
        @test issorted(seq)
        @test all(diff(seq) .> 0)
        for (i, n) in enumerate(m:nmax)
            @test isapprox(seq[i], oblate_cv(m, n, c); rtol=1e-9, atol=1e-10)
        end
    end
end

@testset "Angular functions: public API, parity and ODE residual" begin
    xs = (-0.72, -0.35, 0.2, 0.65)
    cases = ((0, 0, 0.7), (1, 2, 1.1), (2, 4, 1.7))

    for (m, n, c) in cases
        configurations = ((prolate_cv(m, n, c), c^2, x -> prolate_angular_leg(m, n, c, x), x -> prolate_angular_ps(m, n, c, x)), (oblate_cv(m, n, c), -c^2, x -> oblate_angular_leg(m, n, c, x), x -> oblate_angular_ps(m, n, c, x)))

        for (lam, gamma2, leg, ps) in configurations
            for x in xs
                @test all(isfinite, leg(x))
                @test all(isfinite, ps(x))
                leg_residual, leg_d2_error = angular_ode_relative_residual(leg, m, lam, gamma2, x)
                ps_residual, ps_d2_error = angular_ode_relative_residual(ps, m, lam, gamma2, x)
                @test leg_d2_error < 1e-7
                @test ps_d2_error < 1e-7
                @test leg_residual < 5e-8
                @test ps_residual < 5e-8
                @test tuple_isapprox(leg(x), ps(x); rtol=2e-8, atol=1e-10)
            end

            parity = iseven(n - m) ? 1 : -1
            for x in (0.2, 0.35, 0.65, 0.72)
                y_pos, dy_pos = leg(x)
                y_neg, dy_neg = leg(-x)
                @test isapprox(y_neg, parity * y_pos; rtol=1e-10, atol=1e-12)
                @test isapprox(dy_neg, -parity * dy_pos; rtol=1e-10, atol=1e-12)
            end
        end
    end
end

@testset "Known angular precision targets" begin
    targets = ((3, 12, 20.0, 0.2), (8, 16, 40.0, 0.2), (8, 16, 40.0, 0.7))

    for (m, n, c, x) in targets
        residual, d2_error = angular_ode_relative_residual(x -> oblate_angular_leg(m, n, c, x), m, oblate_cv(m, n, c), -c^2, x)
        @test residual < 1e-8
        @test d2_error < 1e-7
    end
end

@testset "Returned first derivatives are coherent" begin
    xs = (-0.72, -0.35, 0.2, 0.65)
    xis = (1.35, 2.2, 3.5)

    for (m, n, c) in ((0, 0, 0.7), (1, 2, 1.1), (2, 4, 1.7))
        angular_functions = (x -> prolate_angular_leg(m, n, c, x), x -> prolate_angular_ps(m, n, c, x), x -> oblate_angular_leg(m, n, c, x), x -> oblate_angular_ps(m, n, c, x))

        for f in angular_functions, x in xs
            @test derivative_relative_error(f, x) < 2e-8
        end

        prolate_radial_functions = (xi -> prolate_radial1(m, n, c, xi), xi -> prolate_radial2(m, n, c, xi))

        for f in prolate_radial_functions, xi in xis
            @test derivative_relative_error(f, xi) < 5e-8
        end

        oblate_radial_functions = (xi -> oblate_radial1(m, n, c, xi), xi -> oblate_radial2(m, n, c, xi))

        for f in oblate_radial_functions, xi in xis
            @test derivative_relative_error(f, xi) < 5e-8
        end
    end
end

@testset "Known derivative precision targets" begin
    targets = ((2, 6, 5.0, 0.7, :oblate_angular), (3, 12, 20.0, -0.6, :prolate_angular), (8, 16, 40.0, 0.2, :oblate_angular), (1, 3, 2.0, 1.35, :oblate_radial2))

    for (m, n, c, x, kind) in targets
        f = kind === :prolate_angular ? (t -> prolate_angular_leg(m, n, c, t)) : kind === :oblate_angular ? (t -> oblate_angular_leg(m, n, c, t)) : kind === :prolate_radial1 ? (t -> prolate_radial1(m, n, c, t)) : kind === :prolate_radial2 ? (t -> prolate_radial2(m, n, c, t)) : kind === :oblate_radial1 ? (t -> oblate_radial1(m, n, c, t)) : (t -> oblate_radial2(m, n, c, t))
        @test derivative_relative_error(f, x) < 1e-10
    end
end

@testset "Expansion coefficients satisfy their recurrence" begin
    for (m, n, c) in ((0, 2, 0.9), (1, 3, 1.4), (2, 5, 2.1))
        pro_lam = prolate_cv(m, n, c)
        pro_dr = SWF.compute_dr2_mix(m, n, c, pro_lam)
        @test coefficient_recurrence_relative_residual(m, n, c, pro_lam, pro_dr) < 1e-12

        obl_lam = oblate_cv(m, n, c)
        obl_dr = SWF.compute_dr2_mix(m, n, im * c, obl_lam)
        @test coefficient_recurrence_relative_residual(m, n, im * c, obl_lam, obl_dr) < 1e-12
    end

    for (m, n, c, max_terms) in ((0, 0, 0.7, 80), (0, 2, 1.3, 100), (2, 5, 4.0, 120), (8, 16, 40.0, 180))
        lam = prolate_cv(m, n, c)
        dr = SWF.compute_dr2_mix(m, n, c, lam, max_terms)
        @test all(isfinite, dr)
        @test coefficient_recurrence_relative_residual(m, n, c, lam, dr) < 1e-12
    end
end

@testset "Known coefficient precision targets" begin
    active_targets = ((2, 6, 5.0, :prolate, 1e-12), (0, 10, 10.0, :prolate, 1e-12), (3, 12, 20.0, :oblate, 1e-12), (8, 16, 40.0, :oblate, 1e-12), (0, 10, 10.0, :oblate, 1e-13))

    for (m, n, c, kind, target_rtol) in active_targets
        lam = kind === :prolate ? prolate_cv(m, n, c) : oblate_cv(m, n, c)
        coeff_c = kind === :prolate ? c : im * c
        dr = SWF.compute_dr2_mix(m, n, coeff_c, lam)
        @test coefficient_recurrence_relative_residual(m, n, coeff_c, lam, dr) < target_rtol
    end

end

@testset "Evenness and c=0 angular limits" begin
    for (m, n, c) in ((0, 0, 0.4), (1, 3, 1.2), (2, 5, 2.0))
        @test isapprox(prolate_cv(m, n, -c), prolate_cv(m, n, c); rtol=1e-12)
        @test isapprox(oblate_cv(m, n, -c), oblate_cv(m, n, c); rtol=1e-12)

        for x in (-0.6, 0.2, 0.7)
            @test tuple_isapprox(prolate_angular_leg(m, n, -c, x), prolate_angular_leg(m, n, c, x); rtol=1e-10)
            @test tuple_isapprox(oblate_angular_leg(m, n, -c, x), oblate_angular_leg(m, n, c, x); rtol=1e-10)
        end
    end

    @test angular_c0_constant_mode_is_constant()
    @test angular_c0_linear_mode_is_linear()

    for m in 0:3, n in m:(m + 3)
        @test prolate_cv(m, n, 0.0) == oblate_cv(m, n, 0.0)
        @test all(tuple_isapprox(prolate_angular_leg(m, n, 0.0, x), oblate_angular_leg(m, n, 0.0, x); rtol=1e-12, atol=1e-13) for x in (-0.8, -0.3, 0.2, 0.7))
    end
end

@testset "Angular orthogonality" begin
    for c in (0.7, 2.0), m in 0:3
        for n in m:(m + 3), k in (n + 1):(m + 4)
            fp = x -> prolate_angular_leg(m, n, c, x)
            gp = x -> prolate_angular_leg(m, k, c, x)
            fo = x -> oblate_angular_leg(m, n, c, x)
            go = x -> oblate_angular_leg(m, k, c, x)

            @test normalized_inner_product(fp, gp) < 2e-9
            @test normalized_inner_product(fo, go) < 2e-9
        end
    end
end

@testset "Known orthogonality precision targets" begin
    active_targets = ((3, 12, 20.0, 1e-12), (8, 16, 40.0, 1e-12), (8, 24, 40.0, 1e-12), (8, 24, 40.0, 1e-13))

    for (m, n, c, target_rtol) in active_targets
        @test normalized_inner_product(x -> prolate_angular_leg(m, n, c, x), x -> prolate_angular_leg(m, n + 1, c, x); nquad=220) < target_rtol
        @test normalized_inner_product(x -> oblate_angular_leg(m, n, c, x), x -> oblate_angular_leg(m, n + 1, c, x); nquad=220) < target_rtol
    end

end

@testset "Near-boundary evaluation" begin
    for (m, n, c) in ((0, 0, 0.7), (1, 2, 1.1), (2, 4, 1.7))
        for x in (-0.999, -0.99, 0.99, 0.999)
            @test all(isfinite, prolate_angular_leg(m, n, c, x))
            @test all(isfinite, oblate_angular_leg(m, n, c, x))
        end

        for xi in (1.001, 1.01, 1.05)
            @test all(isfinite, prolate_radial1(m, n, c, xi))
        end
    end
end

@testset "Radial prolate functions: domain and ODE residual" begin
    r1_xis = (1.35, 2.2, 3.5)
    r2_xis = (2.2, 3.5)
    for (m, n, c) in ((0, 0, 0.7), (0, 2, 1.3), (1, 3, 2.0))
        lam = prolate_cv(m, n, c)
        dr = SWF.compute_dr2_mix(m, n, c, lam)
        r1 = xi -> prolate_radial1(m, n, c, lam, dr, xi)
        r2 = xi -> prolate_radial2(m, n, c, lam, dr, xi)

        for xi in r1_xis
            @test all(isfinite, r1(xi))
            residual, d2_error = radial_ode_relative_residual(r1, m, lam, c, xi)
            @test d2_error < 1e-7
            @test residual < 5e-8
        end

        for xi in r2_xis
            @test all(isfinite, r2(xi))
            residual, d2_error = radial_ode_relative_residual(r2, m, lam, c, xi)
            @test d2_error < 1e-7
            @test residual < 5e-8
        end
    end
end

@testset "Prolate radial Wronskian: Abel identity" begin
    configurations = (((0, 0, 0.7), range(1.4, 8.0; length=16), 5e-3, 5e-3), ((0, 2, 1.3), range(1.4, 8.0; length=16), 5e-4, 5e-4), ((1, 3, 2.0), range(1.4, 8.0; length=16), 5e-4, 5e-4), ((2, 5, 4.0), range(1.4, 8.0; length=16), 1e-4, 1e-4), ((2, 6, 5.0), range(1.8, 10.0; length=16), 1e-10, 1e-10), ((0, 10, 10.0), range(2.0, 15.0; length=16), 1e-10, 1e-10), ((3, 12, 20.0), range(2.5, 20.0; length=16), 1e-10, 1e-10))

    for (case, xis, const_rtol, exact_rtol) in configurations
        m, n, c = case
        scaled_wronskians = map(xi -> scaled_prolate_radial_wronskian(m, n, c, xi), xis)
        reference = scaled_wronskians[length(xis) ÷ 2]
        for value in scaled_wronskians
            @test isapprox(value, reference; rtol=const_rtol, atol=1e-11)
        end
        @test isapprox(reference, inv(c); rtol=exact_rtol, atol=1e-10)
    end

    bulk_xis = range(2.0, 8.0; length=16)
    for (m, n, c) in ((0, 0, 0.7), (0, 2, 1.3), (1, 3, 2.0), (2, 5, 4.0))
        for xi in bulk_xis
            scaled_wronskian = scaled_prolate_radial_wronskian(m, n, c, xi)
            @test isapprox(scaled_wronskian, inv(c); rtol=1e-9, atol=1e-11)
        end
    end
end

@testset "Oblate radial2 near-origin arctan branch" begin
    for (m, n, c) in ((0, 0, 0.7), (0, 2, 1.3), (0, 10, 10.0), (2, 2, 0.7), (2, 6, 5.0), (4, 8, 10.0))
        lam = oblate_cv(m, n, c)
        dr = SWF.compute_dr2_mix(m, n, im * c, lam)

        for xi in (0.001, 0.01, 0.05, 0.1, 0.2)
            R1, dR1 = oblate_radial1(m, n, c, xi)
            R2, dR2 = oblate_radial2(m, n, c, lam, dr, xi)
            W = (1 + xi^2) * (R1 * dR2 - dR1 * R2)

            @test all(isfinite, (R2, dR2))
            @test isapprox(W, inv(c); rtol=1e-10, atol=1e-12)
            @test tuple_isapprox(oblate_radial2(m, n, c, xi), oblate_radial2(m, n, c, lam, xi); rtol=1e-12, atol=1e-12)
            @test tuple_isapprox(oblate_radial2(m, n, c, lam, xi), oblate_radial2(m, n, c, lam, dr, xi); rtol=1e-12, atol=1e-12)
        end

        for xi in (0.001, 0.01, 0.05, 0.1)
            residual, d2_error = oblate_radial_ode_relative_residual(x -> oblate_radial2(m, n, c, lam, dr, x), m, lam, c, xi)
            @test d2_error < 1e-6
            @test residual < 2e-7
        end
    end
end

@testset "Known radial precision targets" begin
    near_boundary_targets = ((0, 0, 0.7, 1.2, 1e-6), (0, 2, 1.3, 1.2, 1e-6), (1, 3, 2.0, 1.2, 1e-6), (2, 5, 4.0, 1.2, 1e-6), (2, 6, 5.0, 1.2, 1e-6), (0, 0, 0.7, 1.4, 1e-6), (0, 2, 1.3, 1.4, 1e-6), (1, 3, 2.0, 1.4, 1e-6), (2, 5, 4.0, 1.4, 1e-6))
    very_near_boundary_targets = ((0, 0, 0.7, 1.1, 1e-10), (0, 2, 1.3, 1.1, 1e-10))
    improved_very_near_boundary_targets = ((1, 3, 2.0, 1.01, 1e-5), (1, 3, 2.0, 1.05, 3e-4), (1, 3, 2.0, 1.075, 6e-4), (1, 3, 2.0, 1.1, 2e-4), (2, 5, 4.0, 1.01, 1e-6), (2, 5, 4.0, 1.05, 3e-5), (2, 5, 4.0, 1.1, 2e-4), (2, 6, 5.0, 1.1, 5e-4), (3, 6, 4.0, 1.1, 2e-5), (4, 8, 10.0, 1.075, 4e-4))
    difficult_very_near_boundary_targets = ((2, 6, 10.0, 1.01, 1e-4), (2, 6, 10.0, 1.02, 8e-4))
    stress_targets = ((8, 16, 40.0, 2.0, 5e-8), (8, 16, 40.0, 5.0, 5e-8), (8, 16, 40.0, 10.0, 5e-8))

    for (m, n, c, xi, target_rtol) in near_boundary_targets
        scaled_wronskian = scaled_prolate_radial_wronskian(m, n, c, xi)
        @test isapprox(scaled_wronskian, inv(c); rtol=target_rtol, atol=1e-11)
    end

    for (m, n, c, xi, target_rtol) in stress_targets
        scaled_wronskian = scaled_prolate_radial_wronskian(m, n, c, xi)
        @test isapprox(scaled_wronskian, inv(c); rtol=target_rtol, atol=1e-11)
    end

    for (m, n, c, xi, target_rtol) in very_near_boundary_targets
        scaled_wronskian = scaled_prolate_radial_wronskian(m, n, c, xi)
        @test isapprox(scaled_wronskian, inv(c); rtol=target_rtol, atol=1e-11)
    end

    for (m, n, c, xi, target_rtol) in improved_very_near_boundary_targets
        scaled_wronskian = scaled_prolate_radial_wronskian(m, n, c, xi)
        @test isapprox(scaled_wronskian, inv(c); rtol=target_rtol, atol=1e-11)
    end

    for (m, n, c, xi, target_rtol) in difficult_very_near_boundary_targets
        scaled_wronskian = scaled_prolate_radial_wronskian(m, n, c, xi)
        @test isapprox(scaled_wronskian, inv(c); rtol=target_rtol, atol=1e-11)
    end
end

@testset "Precomputed parameter overloads are coherent" begin
    xs = [-0.6, -0.1, 0.4, 0.75]
    xis = [1.4, 2.0, 3.0]

    for (m, n, c) in ((0, 2, 0.9), (1, 3, 1.4), (2, 5, 2.1))
        pro_lam = prolate_cv(m, n, c)
        pro_dr = SWF.compute_dr2_mix(m, n, c, pro_lam)
        pro_c2k = SWF.compute_c2k(m, n, pro_dr)
        obl_lam = oblate_cv(m, n, c)
        obl_dr = SWF.compute_dr2_mix(m, n, im * c, obl_lam)
        obl_c2k = SWF.compute_c2k(m, n, obl_dr)

        for x in xs
            @test tuple_isapprox(prolate_angular_leg(m, n, c, x), prolate_angular_leg(m, n, c, pro_lam, pro_dr, x))
            @test tuple_isapprox(prolate_angular_ps(m, n, c, x), prolate_angular_ps(m, n, c, pro_lam, pro_c2k, x))
            @test tuple_isapprox(oblate_angular_leg(m, n, c, x), oblate_angular_leg(m, n, c, obl_lam, obl_dr, x))
            @test tuple_isapprox(oblate_angular_ps(m, n, c, x), oblate_angular_ps(m, n, c, obl_lam, obl_c2k, x))
        end

        for xi in xis
            @test tuple_isapprox(prolate_radial1(m, n, c, xi), prolate_radial1(m, n, c, pro_lam, xi))
            @test tuple_isapprox(prolate_radial1(m, n, c, pro_lam, xi), prolate_radial1(m, n, c, pro_lam, pro_dr, xi))
            @test tuple_isapprox(prolate_radial2(m, n, c, xi), prolate_radial2(m, n, c, pro_lam, xi))
            @test tuple_isapprox(prolate_radial2(m, n, c, pro_lam, xi), prolate_radial2(m, n, c, pro_lam, pro_dr, xi))
            @test tuple_isapprox(oblate_radial1(m, n, c, xi), oblate_radial1(m, n, c, obl_lam, xi))
            @test tuple_isapprox(oblate_radial1(m, n, c, obl_lam, xi), oblate_radial1(m, n, c, obl_lam, obl_dr, xi))
            @test tuple_isapprox(oblate_radial2(m, n, c, xi), oblate_radial2(m, n, c, obl_lam, xi))
            @test tuple_isapprox(oblate_radial2(m, n, c, obl_lam, xi), oblate_radial2(m, n, c, obl_lam, obl_dr, xi))
        end
    end
end
