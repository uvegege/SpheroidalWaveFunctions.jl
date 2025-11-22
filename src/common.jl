"""
    assoc_legendre_Pm(m, N, x) -> (P, dP)

Compute associated Legendre functions Pₙᵐ(x) and their derivatives for n = m, m+1, ..., N.

Uses upward recurrence relations from DLMF for numerical stability.

# Arguments
- `m::Integer`: order (m ≥ 0)
- `N::Integer`: maximum degree (N ≥ m)
- `x::Number`: argument in [-1, 1]

# Returns
- `P::Vector`: values [Pₘᵐ(x), Pₘ₊₁ᵐ(x), ..., Pₙᵐ(x)]
- `dP::Vector`: derivatives [dPₘᵐ/dx, dPₘ₊₁ᵐ/dx, ..., dPₙᵐ/dx]

# Method
- DLMF 14.10.1: Initial value Pₘᵐ(x) = (-1)ᵐ (2m-1)!! (1-x²)^(m/2)
- DLMF 14.10.2: Pₘ₊₁ᵐ(x) = x(2m+1)Pₘᵐ(x)
- DLMF 14.10.3: Upward recurrence for n ≥ m+1
- DLMF 14.10.5: Derivative formula

# Note
For large m, numerical overflow may occur. Consider using a specialized library.
"""
function assoc_legendre_Pm(m, N, x)
    #@assert N >= m "Need N ≥ m."

    P  = zeros(Float64, N - m + 1)
    dP = zeros(Float64, N - m + 1)

    # -------------------------
    # Handle |x| = 1 separately
    # -------------------------
    if abs(x) == 1
        sign = (x == 1 ? +1.0 : -1.0)
        for n = m:N
            idx = n - m + 1
            quotient = 1.0
            if n > m
                for k = (n - m + 1):(n + m)
                    quotient *= k
                end
            end
            P[idx] = sign^(n+m) * quotient 
            dP[idx] = 0.0
        end
        return P, dP
    end

    # -------------------------
    # Step 1: P_m^m(x)
    # DLMF 14.10.1
    # P_m^m(x) = (-1)^m (2m-1)!! (1 - x^2)^(m/2)
    # -------------------------
    if m == 0
        P[1] = 1.0
    else # Not safe for big m.
        # Double factorial (2m - 1)!! 
        df = prod(1:2:(2*m - 1); init = (1 - x^2)^(m/2)) 
        P[1] = (-1)^m * df #* (1 - x^2)^(m/2)
    end

    # Only P_m^m was requested
    if N == m
        # Derivative (using DLMF 14.10.5)
        dP[1] = (m * x * P[1]) / (x^2 - 1)
        return P, dP
    end

    # -------------------------
    # Step 2: P_{m+1}^m(x)
    # DLMF 14.10.2
    # -------------------------
    P[2] = x * (2m + 1) * P[1]

    # -------------------------
    # Step 3: upward recurrence for n ≥ m+1
    # DLMF 14.10.3
    # P_{n+1}^m = ((2n+1)x P_n^m - (n+m)P_{n-1}^m) / (n-m+1)
    # -------------------------
    for n = m+1:N-1
        i = n - m + 1
        P[i+1] = ((2n + 1) * x * P[i] - (n + m) * P[i-1]) / (n - m + 1)
    end

    # -------------------------
    # Step 4: Derivatives
    # DLMF 14.10.5
    # -------------------------
    for n = m:N
        i = n - m + 1
        if n == m
            # Here P_{m-1}^m = 0
            dP[i] = (n * x * P[i]) / (x^2 - 1)
        else
            dP[i] = (n * x * P[i] - (n + m) * P[i-1]) / (x^2 - 1)
        end
    end

    return P, dP
end



"""
    expr_1(n) -> Float64

Compute (2n)!/n! = (n+1)(n+2)⋯(2n).

# Arguments
- `n::Integer`: input value

# Returns
- Value of (2n)!/n!
"""
function expr_1(n)
    a = 1.0
    for i in n+1:2*n
        a *= i
    end
    return a
end

"""
    expr_2(n) -> Float64

Compute (2n+2)!/(2(n+1)!).

# Arguments
- `n::Integer`: input value

# Returns
- Value of (2n+2)!/(2(n+1)!)
"""
function expr_2(n)
    return (2*n+2)*(2*n+1)/(2*(n+1)) * expr_1(n)
end

"""
    expr_3(m, n) -> Float64

Compute the normalization factor for even (n-m) case.

Calculates (n+m)! / ((n-m)/2)! / 2^(n-m).

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number

# Returns
- Normalization factor for even (n-m)
"""
function expr_3(m, n)
    k1 = div(n+m, 2)
    a = 1.0
    for i in k1+1:(n+m)
        a *= i
    end
    for i in 1:div(n-m, 2)
        a /= i
    end
    a /= 2^(n-m)
    return a
end


"""
    expr_4(m, n) -> Float64

Compute the normalization factor for odd (n-m) case.

Calculates (n+m+1)! / ((n-m-1)/2)! / 2^(n-m).

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number

# Returns
- Normalization factor for odd (n-m)
"""
function expr_4(m, n)
    k1 = div(n + m + 1, 2) 
    a = 1.0
    for i in k1+1:(n+m+1)
        a *= i
    end
    for i in 1:div(n-m-1, 2) 
        a /= i
    end
    a /= 2.0^(n - m)
    return a 
end

"""
    αᵣ(m, r, γ²) -> Number

Compute the upper diagonal element of the three-term recurrence relation.

# Arguments
- `m::Integer`: azimuthal order
- `r::Integer`: index
- `γ²::Number`: γ² = c²

# Returns
- Coefficient αᵣ
"""
@inline αᵣ(m, r, γ²) = γ² * ((2 * m + r) + 2) * ((2 * m + r) + 1) / (((2 * m + 2 * r) + 5) * ((2 * m + 2 * r) + 3))

"""
    βᵣ(m, r, γ²) -> Number

Compute the diagonal element of the three-term recurrence relation.

# Arguments
- `m::Integer`: azimuthal order
- `r::Integer`: index
- `γ²::Number`: γ² = c²

# Returns
- Coefficient βᵣ
"""
@inline βᵣ(m, r, γ²) = (m + r) * ((m + r) + 1) + γ² * (2 * (m + r) * ((m + r) + 1) - 2 * m^2 - 1) / (((2 * m + 2 * r) - 1) * ((2 * m + 2 * r) + 3))

"""
    γᵣ(m, r, γ²) -> Number

Compute the lower diagonal element of the three-term recurrence relation.

# Arguments
- `m::Integer`: azimuthal order
- `r::Integer`: index
- `γ²::Number`: γ² = c²

# Returns
- Coefficient γᵣ
"""
@inline γᵣ(m, r, γ²) = γ² * (r * (r - 1)) / (((2 * m + 2 * r) - 3) * ((2 * m + 2 * r) - 1))

"""
    γᵐ(m, r, γ²) -> Number

Alias for βᵣ used in the continued fraction formulation.
"""
@inline γᵐ(m, r, γ²) = βᵣ(m, r, γ²)

"""
    βᵐ(m, r, γ²) -> Number

Compute the product γᵣ(m, r, γ²) * αᵣ(m, r-2, γ²) used in the continued fraction.
"""
@inline βᵐ(m, r, γ²) = γᵣ(m, r, γ²) * αᵣ(m, r - 2, γ²)


"""
    sphericalbesselj_derivative(ν, z) -> Number

Compute the derivative of the spherical Bessel function jᵥ(z) with respect to z.

Uses the recurrence relation: jᵥ'(z) = jᵥ₋₁(z) - (ν+1)/z jᵥ(z).

# Arguments
- `ν::Number`: order
- `z::Number`: argument

# Returns
- Derivative djᵥ/dz
"""
@inline sphericalbesselj_derivative(ν, z) = sphericalbesselj(ν-1, z) - (ν+1)/z * sphericalbesselj(ν, z)

"""
    sphericalbessely_derivative(ν, z) -> Number

Compute the derivative of the spherical Neumann function yᵥ(z) with respect to z.

Uses the recurrence relation: yᵥ'(z) = yᵥ₋₁(z) - (ν+1)/z yᵥ(z).

# Arguments
- `ν::Number`: order
- `z::Number`: argument

# Returns
- Derivative dyᵥ/dz
"""
@inline sphericalbessely_derivative(ν, z) = sphericalbessely(ν-1, z) - (ν+1)/z * sphericalbessely(ν, z)




"""
    Nmn(m, n, dr) -> Float64

Compute the normalization constant Nₘₙ(c) for spheroidal wave functions.

Calculates the normalization constant using the expansion coefficients and
factorial recursion relations.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `dr::Vector`: expansion coefficients

# Returns
- Normalization constant Nₘₙ(c)

# Note
The normalization ensures proper orthogonality and norm conditions for the
spheroidal wave functions.
"""
function Nmn(m, n, dr)
    
    init = iseven(n - m) ? 0 : 1
    nmax = length(dr)
    r = init 
    fact = 1.0
    for i in (init + 1):(2*m + init)
        fact *= i
    end
    
    sum = dr[1]^2 * fact / (2*m+2*r+1)
    for i in 2:nmax
        r += 2 
        fact *= (2*m+r) * (2*m + r - 1)/(r*(r-1))
        sum += dr[i]^2 * fact / (2*m+2*r+1)
    end
    
    return 2*sum
end


