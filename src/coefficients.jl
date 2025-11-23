"""
    compute_dr2_mix(m, n, c, λ, max_terms=25 + div(n - m, 2) + round(Int, abs(c)))

Compute the expansion coefficients dᵣ for the spheroidal wave function using a mixed method.

Combines backward and forward recursion to compute the expansion coefficients in the
Legendre function expansion. The coefficients satisfy a three-term recurrence relation.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter
- `λ::Number`: characteristic value
- `max_terms::Int`: maximum number of expansion terms (default: adaptive)

# Returns
- Vector of normalized expansion coefficients dᵣ

# Method
Uses backward recursion from the tail and forward recursion from the start, joining
at an optimal point to maintain numerical stability.
"""
function compute_dr2_mix(m, n, c, λ, max_terms = 25 + div(n - m, 2) + round(Int, abs(c)))

    γ² = real(c^2)
    is_even = iseven(n - m)
    point = div((n-m+2), 2) - (floor(Int, log10(abs(c))) + 1)
    #point = div((n-m), 2) - (floor(Int, log10(abs(c))) + 1)

    if point <= 1
        point = 1
    end
    dr = zeros(max_terms)
    r_values = is_even ? (0:2:2*(max_terms-1)) : (1:2:2*max_terms-1)
    
    start_idx = is_even ? length(r_values) : length(r_values) - 1
    dr[start_idx] = 1.0
    N = 1.0
    # Compute coefficients downwards using continued fraction
    for i in start_idx-1:-1:max(point, 1)
        r = r_values[i+1]  # Current r
        r_prev = r_values[i]  # r-2
        if i == start_idx-1
            N = 1.0
        else
            β_m = βᵐ(m, r, γ²)
            γ_m = γᵐ(m, r, γ²)
            N = β_m / (γ_m - λ - N)
        end
        α_val = αᵣ(m, r_prev, γ²) 
        dr[i] = -(α_val / N) * dr[i+1] # Compute d_{r-2} from d_r
    end

    # Compute coefficients with forward marching
    if point > 1
        dr[1] = 1.0
        f = 1.0
        r = r_values[1]
        value = -((βᵣ(m, r, γ²) - λ) * dr[1])/αᵣ(m, r, γ²)
        f = value
        if point > 2
            r = r_values[1]
            dr[2] = value
            f = dr[2]
            for i in 2:point-1
                r = r_values[i]
                value = -((βᵣ(m, r, γ²) - λ) * dr[i] + γᵣ(m, r, γ²)*dr[i-1])/αᵣ(m, r, γ²)
                f = value
                if i != point-1
                    dr[i+1] = value
                end
            end
        end

        # Normalize the values
        C = dr[point] / f
        for i in 1:point-1
            dr[i] *= C
        end
    end

    x = compute_normalization_sum(dr, r_values, m, is_even)
    s = compute_normalization_constant(m, n, is_even) / x

    dr .*= s
    
    return dr
end


"""
    compute_normalization_sum(dr, r_values, m, is_even)

Compute the normalization sum for the expansion coefficients.

# Arguments
- `dr::Vector`: expansion coefficients
- `r_values::Range`: range of r indices
- `m::Integer`: azimuthal order
- `is_even::Bool`: true if (n-m) is even

# Returns
- Normalization sum value
"""
function compute_normalization_sum(dr, r_values, m, is_even)
    x = 0.0
    a = 1.0
    for (i, r) in enumerate(r_values)
        if is_even
            if r > 0
                a *= (-(2*m + r - 1) * (2*m + r)) / (4 * (m + r/2) * (r/2))
            else
                a = expr_1(m) # factorial(2m) / factorial(m)
            end
        else
            if r > 1
                a *= (-(2*m + r) * (2*m + r + 1)) / (4 * (m + (r+1)/2) * ((r-1)/2))
            else
                a = expr_2(m) # factorial(2m + 2) / (2 * factorial(m + 1))
            end
        end
        x += dr[i] * a
    end
    return x
end

"""
    compute_normalization_constant(m, n, is_even)

Compute the normalization constant based on the parity of (n-m).

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `is_even::Bool`: true if (n-m) is even

# Returns
- Normalization constant
"""
function compute_normalization_constant(m, n, is_even)
    if is_even
        return (-1)^((n-m) ÷ 2) * expr_3(m, n)
    else
        return (-1)^((n-m-1) ÷ 2) * expr_4(m, n)
    end
end


"""
    pochhammer(a, k)

Compute the Pochhammer symbol (rising factorial) (a)ₖ = a(a+1)(a+2)⋯(a+k-1).

# Arguments
- `a::Number`: base value
- `k::Integer`: number of terms

# Returns
- Pochhammer symbol (a)ₖ
"""
function pochhammer(a, k)
    result = 1.0
    for i in 0:k-1
        result *= (a + i)
    end
    return result
end

"""
    compute_c2k(m, n, dr)

Compute the power series expansion coefficients c₂ₖ from the Legendre expansion coefficients dᵣ.

Transforms the expansion coefficients from the Legendre function basis to the power series
basis (1-x²)ᵏ, which is useful for certain computational approaches.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `dr::Vector`: Legendre expansion coefficients

# Returns
- Vector of power series coefficients c₂ₖ

# Notes
The recurrence relations used are:
- (2m+r)!/r! → (2m+r+1)(2m+r+2)/((r+1)(r+2))
- (-r/2)ₖ → (-r/2 + 1) / (-r/2 + k - 1)
- (m+r/2+0.5)ₖ → (m + r/2 + 0.5 + k)/(m + r/2 + 0.5)
"""
function compute_c2k(m, n, dr)

    N = length(dr)
    c2k = zeros(N)
    is_even = iseven(n - m)
    kfact = N < 80 ? 1.0 : 1e-200 # used on scipy, i had problems when length(dr) is big
    for i in 1:m # 1/((m+k)!/k!)
        kfact *= i 
    end
    for k in 0:N-1
        sum_val = 0.0
        if is_even
            r = 2*k
            fact = 1.0
            for i in r+1:r+2*m # (2m+r)!/r! -> (r+1)*(r+2)*...*(r+2m)
                fact *= i 
            end
            pk = pochhammer(-r/2, k) # Initialize  (-r/2)ₖ 
            pmk = pochhammer(m+r/2+1/2, k) # Initialize (m+r/2+1/2, k)ₖ 
            sum_val += dr[k+1] * fact * pk * pmk
            current_term = 0.0
            for k2 in k+1:N-1
                r = 2*(k2-1) # sum over r = 2k, 2k+2, 2k+4, ...
                fact *= ((2*m+r+1)*(2*m+r+2)) / ((r+1)*(r+2))
                pk *= (-r/2 - 1) / (-r/2 + k - 1) # (-r/2 + 1) ?
                pmk *= (m + r/2 + 0.5 + k)/(m + r/2 + 0.5)
                term = dr[k2+1] * fact * pk * pmk
                sum_val += term
                abs(current_term - sum_val) < abs(sum_val) * 1e-14 && break
                current_term = sum_val
            end
        else
            r = 2*k+1
            fact = 1.0
            for i in r+1:r+2*m # (2m+r)!/r! -> (r+1)*(r+2)*...*(r+2m)
                fact *= i 
            end
            pk = pochhammer(-(r-1)/2, k) # Initialize  (-(r-1)/2)ₖ
            pmk = pochhammer(m+r/2+1, k) # Initialize -(m+r/2+1/2, k)ₖ
            current_term = 0.0
            sum_val += dr[k+1] * fact * pk * pmk
            for k2 in k+1:N-1
                r = 2*(k2-1)+1 # sum over r = 2k+1, 2k+3, 2k+5, ...
                fact *= (2*m+r+1)*(2*m+r+2)/((r+1)*(r+2))
                pk *= (-r/2 - 1/2) / (-r/2 - 0.5 + k)
                pmk *= (m + r/2 + 1 + k)/(m + r/2 + 1)
                term = dr[k2+1] * fact * pk * pmk
                sum_val += term
                abs(current_term - sum_val) < abs(sum_val) * 1e-14 && break
                current_term = sum_val
            end
        end
        c2k[k+1] = sum_val / (2^m) / kfact
        kfact *= (m + k + 1) * (k + 1)
    end
    return c2k
end
