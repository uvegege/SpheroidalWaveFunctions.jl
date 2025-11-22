"""
    Fmn(m, n, dr)

Compute the normalization factor Fₘₙ for radial spheroidal wave functions.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `dr::Vector`: expansion coefficients

# Returns
- Normalization factor Fₘₙ
"""
function Fmn(m, n, dr)
    init = iseven(n-m) ? 0 : 1
    nmax = length(dr)
    sum = 0.0
    fac = 1.0
    for i in init+1:2*m + init
        fac *= i
    end
    sum = dr[1] * fac
    r = init
    for i in 2:nmax
        r += 2
        fac *= (2*m+r) * (2*m + r - 1)/(r*(r-1))
        sum += dr[i] * fac
    end
    return sum
end

"""
    spheroidal_rad_1(m, n, c, dr, ξ)

Compute the radial spheroidal wave function of the first kind using spherical Bessel functions.

The first kind radial function is regular at ξ = 1 and is expressed as a series
of spherical Bessel functions jₙ.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter
- `dr::Vector`: expansion coefficients
- `ξ::Number`: radial argument (ξ ≥ 1)

# Returns
- `R`: value of the radial function Rₘₙ⁽¹⁾(c, ξ)
- `dR`: derivative dRₘₙ⁽¹⁾/dξ
"""
function spheroidal_rad_1(m, n, c, dr, ξ)
    #λ =  find_eigenvalue(m, n, c)
    #dr = compute_dr2_mix(m, n, c, λ)
    Fmn_val = Fmn(m, n, dr)
    factor = 1/Fmn_val * (1 - 1/real(ξ^2))^(m/2)
    dfactor = m * (1 - 1/real(ξ^2)) ^((m-2)/2) * abs(ξ)^(-3)  / Fmn_val
    init = iseven(n-m) ? 0 : 1
    nmax = length(dr)

    fact = 1.0
    for i in init+1:2*m + init
        fact *= i
    end
    r = init
    v = (-1)^((r - (n-m))/2) * dr[1] * fact 
    S = v * sphericalbesselj(m+r, abs(c)*abs(ξ))
    dS = v * sphericalbesselj_derivative(m+r, abs(c)*abs(ξ)) * abs(c)
    for i in 2:nmax
        r += 2
        fact *= (2*m+r) * (2*m + r - 1)/(r*(r-1))
        v = (-1)^((r - (n-m))/2) * dr[i] * fact 
        S += v * sphericalbesselj(m+r, abs(c)*abs(ξ))
        dS += v * sphericalbesselj_derivative(m+r, abs(c)*abs(ξ)) * abs(c)
    end

    dS = factor * dS + dfactor * S
    S = S * factor
    return S, dS
end

"""
    spheroidal_rad_2(m, n, c, dr, ξ)

Compute the radial spheroidal wave function of the second kind using spherical Neumann functions.

The second kind radial function is irregular at ξ = 1 and is expressed as a series
of spherical Neumann functions yₙ.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter
- `dr::Vector`: expansion coefficients
- `ξ::Number`: radial argument (ξ ≥ 1)

# Returns
- `R`: value of the radial function Rₘₙ⁽²⁾(c, ξ)
- `dR`: derivative dRₘₙ⁽²⁾/dξ
"""
function spheroidal_rad_2(m, n, c, dr, ξ)
    #λ =  find_eigenvalue(m, n, c)
    #dr = compute_dr2_mix(m, n, c, λ)
    Fmn_val = Fmn(m, n, dr)
    factor = 1/Fmn_val * (1 - 1/real(ξ^2))^(m/2)
    dfactor = m * (1 - 1/real(ξ^2)) ^((m-2)/2) * abs(ξ)^(-3)  / Fmn_val
    init = iseven(n-m) ? 0 : 1
    nmax = length(dr)

    fact = 1.0
    for i in init+1:2*m + init
        fact *= i
    end
    r = init
    v = (-1)^((r - (n-m))/2) * dr[1] * fact 
    S = v * sphericalbessely(m+r, abs(c)*abs(ξ))
    dS = v * sphericalbessely_derivative(m+r, abs(c)*abs(ξ)) * abs(c)
    for i in 2:nmax
        r += 2
        fact *= (2*m+r) * (2*m + r - 1)/(r*(r-1))
        v = (-1)^((r - (n-m))/2) * dr[i] * fact 
        S += v * sphericalbessely(m+r, abs(c)*abs(ξ))
        dS += v * sphericalbessely_derivative(m+r, abs(c)*abs(ξ)) * abs(c)
    end

    dS = factor * dS + dfactor * S
    S = S * factor
    return S, dS
end

"""
    kmn1(m, n, c)

Compute the normalization constant kₘₙ for the first radial function.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter

# Returns
- Normalization constant kₘₙ

# Reference
DLMF: (30.11.10) / (30.11.11)
"""
function kmn1(m, n, c)
    λ =  find_eigenvalue(m, n, c)
    dr = compute_dr2_mix(m, n, c, λ)
    v = 1.0
    if iseven(n-m)
        #factorial(m+n)/ (factorial(div(n-m, 2)) * factorial(div(n+m, 2)))
        k1 = div(n+m, 2)
        v = 1.0
        for i in k1+1:(n+m)
            v *= i
        end
        for i in 1:div(n-m, 2)
            v /= i
        end 
        v = v * (2*m+1) * Fmn(m, n, dr) / (2^(m+n)) / dr[1] / c^m
        for i in 1:m # 1/factorial(m)
            v /= i
        end

    else
        k1 = div(n + m + 1, 2) 
        v = 1.0
        for i in k1+1:(n+m+1)
            v *= i
        end
        for i in 1:div(n-m-1, 2) 
            v /= i
        end
        v = v * (2*m+3) * Fmn(m, n, dr) / (2^(m+n)) / dr[1] / c^(m+1)
        for i in 1:m # 1/factorial(m)
            v /= i
        end
    end
    return v
end

#TODO: not ok when ξ > 3
function spheroidal_rad_1_v1(m, n, c, ξ)
    kfactor = 1/kmn1(m, n, c)
    P, dP = spheroidal_ang_1.(m, n, c, ξ)
    return -kfactor * P, -kfactor * dP
end

#TODO: not ok when ξ > 3
function spheroidal_rad_1_v2(m, n, c, ξ)
    kfactor = 1/kmn1(m, n, c)
    P, dP = spheroidal_ang_2.(m, n, c, ξ)
    return -kfactor * P, -kfactor * dP
end
