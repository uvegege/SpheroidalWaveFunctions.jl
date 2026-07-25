"""
Special radial-function evaluation paths.

These functions implement auxiliary expansions for the radial spheroidal
function of the second kind in parameter regions where the direct
spherical-Neumann series is poorly conditioned.  The formulas follow the
methods described by Adelman, Gumerov, and Duraiswami in "Software for
Computing the Spheroidal Wave Functions Using Arbitrary Precision
Arithmetic", especially the prolate Legendre-Q representation and the
oblate arctan power-series representation near the singular/small-argument
regions.
"""

function legendre_pqmns(m, N, x)
    x <= 1 && throw(DomainError(x, "Legendre Q requires x > 1"))

    root = sqrt(x^2 - 1)
    P = zeros(m + 2, N + 2)
    Q = zeros(m + 2, N + 2)

    P[1, 1] = 1.0
    Q[1, 1] = 0.5 * log((x + 1) / (x - 1))
    P[1, 2] = x
    Q[1, 2] = x * Q[1, 1] - 1

    for n in 1:N
        P[1, n+2] = ((2n + 1) * x * P[1, n+1] - n * P[1, n]) / (n + 1)
        Q[1, n+2] = ((2n + 1) * x * Q[1, n+1] - n * Q[1, n]) / (n + 1)
    end

    for order in 0:m, n in order+1:N+1
        P[order+2, n+1] = ((n - order) * x * P[order+1, n+1] - (n + order) * P[order+1, n]) / root
        Q[order+2, n+1] = ((n - order) * x * Q[order+1, n+1] - (n + order) * Q[order+1, n]) / root
    end

    p = P[m+1, 1:N+1]
    q = Q[m+1, 1:N+1]
    dp = [(-(n + 1) * x * P[m+1, n+1] + (n - m + 1) * P[m+1, n+2]) / (x^2 - 1) for n in 0:N]
    dq = [(-(n + 1) * x * Q[m+1, n+1] + (n - m + 1) * Q[m+1, n+2]) / (x^2 - 1) for n in 0:N]

    if m > 0
        Qm = zeros(m + 1, m + 1)
        dQm = zeros(m + 1, m + 1)
        Qm[1, 1] = Q[1, 1]
        dQm[1, 1] = -1 / (x^2 - 1)
        Qm[1, 2] = Q[1, 2]
        dQm[1, 2] = Qm[1, 1] + x * dQm[1, 1]

        for n in 1:m-1
            Qm[1, n+2] = ((2n + 1) * x * Qm[1, n+1] - n * Qm[1, n]) / (n + 1)
            dQm[1, n+2] = ((2n + 1) * (Qm[1, n+1] + x * dQm[1, n+1]) - n * dQm[1, n]) / (n + 1)
        end

        for order in 0:m-1, n in 0:m
            q2 = ((n * (n + 1) + order^2 / (x^2 - 1)) * Qm[order+1, n+1] - 2*x * dQm[order+1, n+1]) / (x^2 - 1)
            Qm[order+2, n+1] = root * dQm[order+1, n+1] - order * x * Qm[order+1, n+1] / root
            dQm[order+2, n+1] = root * q2 + (1 - order) * x * dQm[order+1, n+1] / root + order * Qm[order+1, n+1] / root^3
        end

        for n in 0:m-1
            q[n+1] = Qm[m+1, n+1]
            dq[n+1] = dQm[m+1, n+1]
        end
    end

    return p, dp, q, dq
end

function compute_dr_negative(m, n, c, λ, dr)
    is_even = iseven(n - m)
    stop_r = is_even ? -2*m : -2*m + 1
    r0 = is_even ? 0 : 1
    nterms = div(r0 - stop_r, 2)
    nterms == 0 && return zeros(0)

    γ² = c^2
    A = zeros(nterms, nterms)
    b = zeros(nterms)

    for i in 1:nterms
        r = stop_r + 2 * (i - 1)
        i > 1 && (A[i, i-1] = γᵣ(m, r, γ²))
        A[i, i] = βᵣ(m, r, γ²) - λ
        i < nterms ? (A[i, i+1] = αᵣ(m, r, γ²)) : (b[i] = -αᵣ(m, r, γ²) * dr[1])
    end

    return A \ b
end

function negative_regularized_denominator(m, c, λ, r, min_r)
    γ² = c^2
    first_r = iseven(r - min_r) ? min_r : min_r - 1
    value = βᵣ(m, first_r, γ²) - λ

    for rr in first_r+2:2:r-2
        value = βᵣ(m, rr, γ²) - λ - γᵣ(m, rr, γ²) * αᵣ(m, rr - 2, γ²) / value
    end

    return βᵣ(m, r, γ²) - λ - γᵣ(m, r, γ²) * αᵣ(m, r - 2, γ²) / value
end

function compute_dr_negative_regularized(m, n, c, λ, dr, dr_negative; min_extra=220, rtol=1e-14)
    is_even = iseven(n - m)
    stop_r = is_even ? -2*m : -2*m + 1
    start_r = stop_r - 2
    min_r = start_r - 2 * min_extra
    coeffs = Float64[]
    anchor = m == 0 ? dr[1] : dr_negative[1]
    factor = is_even ? c^2 / ((2*m - 1) * (2*m + 1)) : -c^2 / ((2*m - 1) * (2*m - 3))

    push!(coeffs, anchor * factor / negative_regularized_denominator(m, c, λ, start_r, min_r))

    for r in start_r-2:-2:min_r
        value = coeffs[end] * (-αᵣ(m, r, c^2) / negative_regularized_denominator(m, c, λ, r, min_r))
        push!(coeffs, value)
        abs(value) < rtol * maximum(abs, coeffs) && break
    end

    return coeffs
end

function product_to(n)
    value = 1.0
    for i in 1:n
        value *= i
    end
    return value
end

function kmn2_prolate(m, n, c, dr, dr_negative)
    F = Fmn(m, n, dr)

    if iseven(n - m)
        dstop = m == 0 ? dr[1] : dr_negative[1]
        return product_to(div(n - m, 2)) * product_to(div(m + n, 2)) * dstop * F * 2.0^(n - m) * product_to(2*m) / ((2*m - 1) * product_to(m) * product_to(m + n) * c^(m - 1))
    else
        dstop = dr_negative[1]
        return -product_to(div(n - m - 1, 2)) * product_to(div(m + n + 1, 2)) * dstop * F * 2.0^(n - m) * product_to(2*m) / ((2*m - 3) * (2*m - 1) * product_to(m) * product_to(m + n + 1) * c^(m - 2))
    end
end

q_index(l) = l >= 0 ? l + 1 : -l

function prolate_radial2_legendreQ(m, n, c, λ, dr, ξ; rtol=1e-11, window=3)
    dr_negative = compute_dr_negative(m, n, c, λ, dr)
    dr_regularized = compute_dr_negative_regularized(m, n, c, λ, dr, dr_negative)
    return prolate_radial2_legendreQ(m, n, c, dr, dr_negative, dr_regularized, ξ; rtol=rtol, window=window)
end

function prolate_radial2_legendreQ(m, n, c, dr, dr_negative, dr_regularized, ξ; rtol=1e-11, window=3)
    is_even = iseven(n - m)
    stop_r = is_even ? -2*m : -2*m + 1
    r0 = is_even ? 0 : 1
    max_degree = max(m + r0 + 2 * (length(dr) - 1), m - (stop_r - 2 * length(dr_regularized)) - 1, m + 8)
    pm, pd, qm, qd = legendre_pqmns(m, max_degree, ξ)

    S = 0.0
    dS = 0.0
    terms = zeros(length(dr_negative) + length(dr) + length(dr_regularized))
    dterms = similar(terms)
    i = 0

    for (j, r) in enumerate(stop_r:2:r0-2)
        i += 1
        terms[i] = dr_negative[j] * qm[q_index(m + r)]
        dterms[i] = dr_negative[j] * qd[q_index(m + r)]
        S += terms[i]
        dS += dterms[i]
    end

    r = r0
    for j in eachindex(dr)
        i += 1
        terms[i] = dr[j] * qm[m+r+1]
        dterms[i] = dr[j] * qd[m+r+1]
        S += terms[i]
        dS += dterms[i]
        r += 2
    end

    r = stop_r - 2
    for j in eachindex(dr_regularized)
        degree = -r - m - 1
        i += 1
        terms[i] = dr_regularized[j] * pm[degree+1]
        dterms[i] = dr_regularized[j] * pd[degree+1]
        S += terms[i]
        dS += dterms[i]
        r -= 2
    end

    tail_start = max(1, length(terms) - window + 1)
    converged = maximum(abs, @view terms[tail_start:end]) / (abs(S) + eps()) < rtol && maximum(abs, @view dterms[tail_start:end]) / (abs(dS) + eps()) < rtol
    k2 = kmn2_prolate(m, n, c, dr, dr_negative)

    return S / k2, dS / k2, converged
end

function use_prolate_radial2_legendreQ(m, n, c, ξ)
    δ = ξ - 1
    ρ = c / (n + 1)
    within(limit) = δ <= limit + 16*eps(float(ξ))
    m == 0 && return iseven(n) && within(0.15)
    m == 1 && return ρ <= 1.5 && within(0.075)
    m == 2 && return ρ <= 1.0 ? within(0.10) : ρ <= 1.5 && within(0.02)
    m == 3 && return ρ <= 1.5 && within(0.10)
    return m == 4 && ρ <= 1.5 && within(0.10)
end

function compute_qstar_alpha(c2k; nterms=length(c2k))
    B = zeros(nterms)
    for n in 0:nterms-1
        for k in 0:n
            n - k + 1 <= nterms && (B[n+1] += c2k[k+1] * c2k[n-k+1])
        end
    end

    A = zeros(nterms)
    A[1] = inv(B[1])
    for n in 1:nterms-1
        for k in 0:n-1
            A[n+1] -= A[k+1] * B[n-k+1] / B[1]
        end
    end

    alpha = zeros(nterms)
    factor = 1.0
    for r in 0:nterms-1
        r > 0 && (factor *= r)
        alpha[r+1] = A[r+1] * factor
    end

    return alpha
end

function compute_oblate_qstar_even(m, c, k1, c2k)
    alpha = compute_qstar_alpha(c2k)
    value = 0.0 + 0.0im

    for r in 0:m
        value += alpha[r+1] * product_to(2*m - 2r) / (product_to(r) * (2.0^(m - r) * product_to(m - r))^2)
    end

    return ((1.0im)^(-m) * k1)^2 * value / c
end

function compute_oblate_h2r_even(m, c, qstar, k1, c2k, r)
    value = 0.0 + 0.0im

    for k in max(0, r - m + 1):length(c2k)-1
        value += c2k[k+1] * (m + 2k) * product_to(m + k - 1) / (product_to(m + k - 1 - r) * product_to(r))
    end

    return -2 * qstar * value / ((1.0im)^(-m) * k1)
end

function compute_oblate_B2r_even(m, n, c, λ, qstar, k1, c2k, R1_0; nterms=length(c2k))
    B = zeros(ComplexF64, nterms)
    B[1] = c / R1_0 - qstar * R1_0

    for r in 0:nterms-2
        α2r = (2r + 2) * (2r + 3)
        β2r = (2r + 1) * (2r - 2*m + 2) + m * (m - 1) - λ
        h2r = compute_oblate_h2r_even(m, c, qstar, k1, c2k, r)
        B[r+2] = (h2r - β2r * B[r+1] - (r >= 1 ? c^2 * B[r] : 0.0)) / α2r
    end

    return B
end

function oblate_radial1_value_at_zero(m, n, c, dr)
    value = first(spheroidal_rad_1(m, n, im * c, dr, 0.0 + 0.0im))
    isfinite(real(value)) && return value

    for ξ in (1e-10, 1e-8, 1e-6, 1e-4)
        value = first(spheroidal_rad_1(m, n, im * c, dr, im * ξ))
        isfinite(real(value)) && return value
    end

    return value
end

function eval_oblate_g_even(m, B, ξ; rtol=1e-14)
    ξ2 = ξ^2
    power = 1.0
    poly = 0.0 + 0.0im
    dpoly = 0.0 + 0.0im

    for r in eachindex(B)
        k = r - 1
        term = B[r] * power
        poly += term
        k > 0 && (dpoly += 2k * B[r] * ξ^(2k - 1))
        abs(term) < rtol * abs(poly) && break
        power *= ξ2
    end

    factor = (1 + ξ2)^(-m / 2)
    dfactor = -m * ξ * (1 + ξ2)^(-m / 2 - 1)
    g = ξ * factor * poly
    dg = factor * poly + ξ * dfactor * poly + ξ * factor * dpoly

    return g, dg
end

function oblate_radial2_arctan_even_data(m, n, c, λ, dr)
    c2k = compute_c2k(m, n, dr)
    k1 = kmn1(m, n, im * c)
    qstar = compute_oblate_qstar_even(m, c, k1, c2k)
    R1_0 = oblate_radial1_value_at_zero(m, n, c, dr)
    B = compute_oblate_B2r_even(m, n, c, λ, qstar, k1, c2k, R1_0)
    return qstar, B
end

function oblate_radial2_arctan_even(m, n, c, dr, qstar, B, ξ)
    R1, dR1 = spheroidal_rad_1(m, n, im * c, dr, im * ξ)
    g, dg = eval_oblate_g_even(m, B, ξ)
    h0 = atan(ξ) - π / 2
    R = qstar * R1 * h0 + g
    dR = qstar * dR1 * h0 + qstar * R1 / (1 + ξ^2) + dg

    return real(R / c^2), real(dR / c^2)
end

function oblate_radial2_arctan_even(m, n, c, λ, dr, ξ)
    qstar, B = oblate_radial2_arctan_even_data(m, n, c, λ, dr)
    return oblate_radial2_arctan_even(m, n, c, dr, qstar, B, ξ)
end

function use_oblate_radial2_arctan(m, n, c, ξ)
    return iseven(m) && iseven(n - m) && c / (n + 1) <= 1.5 && ξ <= 0.25
end

function spheroidal_rad_2(m, n, c, λ, dr, ξ)
    if !isreal(c) && !isreal(ξ)
        c_oblate = abs(im * c)
        ξ_oblate = abs(im * ξ)
        if use_oblate_radial2_arctan(m, n, c_oblate, ξ_oblate)
            R, dR = oblate_radial2_arctan_even(m, n, c_oblate, λ, dr, ξ_oblate)
            isfinite(R) && isfinite(dR) && return R, dR
        end
        return spheroidal_rad_2(m, n, c, dr, ξ)
    end

    use_prolate_radial2_legendreQ(m, n, c, ξ) || return spheroidal_rad_2(m, n, c, dr, ξ)
    R, dR, converged = prolate_radial2_legendreQ(m, n, c, λ, dr, ξ)
    return isfinite(R) && isfinite(dR) ? (R, dR) : spheroidal_rad_2(m, n, c, dr, ξ)
end

function spheroidal_rad_2(m, n, c, λ, dr, ξs::AbstractArray)
    if !isreal(c)
        values = Vector{Tuple{Float64, Float64}}(undef, length(ξs))
        c_oblate = abs(im * c)
        ξs_oblate = abs.(im .* ξs)
        use_special = any(ξ -> use_oblate_radial2_arctan(m, n, c_oblate, ξ), ξs_oblate)
        qstar, B = use_special ? oblate_radial2_arctan_even_data(m, n, c_oblate, λ, dr) : (0.0 + 0.0im, ComplexF64[])

        for (i, ξ) in enumerate(ξs)
            ξ_oblate = abs(im * ξ)
            if use_oblate_radial2_arctan(m, n, c_oblate, ξ_oblate)
                R, dR = oblate_radial2_arctan_even(m, n, c_oblate, dr, qstar, B, ξ_oblate)
                values[i] = isfinite(R) && isfinite(dR) ? (R, dR) : spheroidal_rad_2(m, n, c, dr, ξ)
            else
                values[i] = spheroidal_rad_2(m, n, c, dr, ξ)
            end
        end
        return reshape(values, axes(ξs))
    end

    values = Vector{Tuple{Float64, Float64}}(undef, length(ξs))
    use_special = any(ξ -> use_prolate_radial2_legendreQ(m, n, c, ξ), ξs)
    dr_negative = use_special ? compute_dr_negative(m, n, c, λ, dr) : zeros(0)
    dr_regularized = use_special ? compute_dr_negative_regularized(m, n, c, λ, dr, dr_negative) : zeros(0)

    for (i, ξ) in enumerate(ξs)
        if use_prolate_radial2_legendreQ(m, n, c, ξ)
            R, dR, converged = prolate_radial2_legendreQ(m, n, c, dr, dr_negative, dr_regularized, ξ)
            values[i] = isfinite(R) && isfinite(dR) ? (R, dR) : spheroidal_rad_2(m, n, c, dr, ξ)
        else
            values[i] = spheroidal_rad_2(m, n, c, dr, ξ)
        end
    end

    return reshape(values, axes(ξs))
end
