"""
    spheroidal_ang_1(m, n, dr, x)

Compute the angular spheroidal wave function using Legendre function expansion (method 1).

Evaluates the spheroidal function as:
- Sₘₙ(x, γ²) = ∑ₖ d₂ₖ(c) Pₘ₊₂ₖᵐ(x)  for even (n-m)
- Sₘₙ(x, γ²) = ∑ₖ d₂ₖ₊₁(c) Pₘ₊₂ₖ₊₁ᵐ(x)  for odd (n-m)

where Pₙᵐ are associated Legendre functions.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `dr::Vector`: expansion coefficients dᵣ
- `x::Number`: argument in [-1, 1]

# Returns
- `S`: value of the angular function Sₘₙ(x)
- `dS`: derivative dSₘₙ/dx
"""
function spheroidal_ang_1(m, n, dr, x)
    #λ =  find_eigenvalue(m, n, c)
    #dr = compute_dr2_mix(m, n, c, λ)
    v = iseven(n-m) ? 0 : 1
    S = 0.0
    dS = 0.0
    plms = assoc_legendre_Pm(m, m+2*(length(dr)-1)+v, x)

    for k in eachindex(dr)
        coef = dr[k]
        Plm  = plms[1][2*k-1+v]
        dPlm = plms[2][2*k-1+v]
        S += coef * Plm
        dS += coef * dPlm
    end
    s = iseven(m) ? 1 : -1
    return s*S, s*dS
end


"""
    spheroidal_ang_2(m, n, c2k, x)

Compute the angular spheroidal wave function using power series expansion (method 2).

Evaluates the spheroidal function as:
- Sₘₙ(x, γ²) = (-1)ᵐ (1 - x²)^(m/2) ∑ₖ c₂ₖ (1 - x²)ᵏ  for even (n-m)
- Sₘₙ(x, γ²) = (-1)ᵐ x(1 - x²)^(m/2) ∑ₖ c₂ₖ (1 - x²)ᵏ  for odd (n-m)

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c2k::Vector`: power series coefficients c₂ₖ
- `x::Number`: argument in [-1, 1]

# Returns
- `S`: value of the angular function Sₘₙ(x)
- `dS`: derivative dSₘₙ/dx
"""
function spheroidal_ang_2(m, n, c2k, x)
    #λ =  find_eigenvalue(m, n, c)
    #dr = compute_dr2_mix(m, n, c, λ)
    #c2k = compute_c2k(m, n, c, λ, dr)
    S = 0.0
    dS = 0.0
    for r in eachindex(c2k)
        k = r-1
        coef = c2k[r]
        S += coef* (1-x^2)^(k)
        dS += coef * (-2*k*x) * (1-x^2)^(k-1)
    end
    if iseven(n-m)
        dfactor = -m*x*(1-x^2)^(m/2-1)
        factor = (1-x^2)^(m/2)
        dS = (-1)^m * (factor * dS + dfactor * S)
        S = (-1)^m * factor * S
    else
        dfactor = (1-x^2)^(m/2 - 1) * (1 - (m+1) * x^2)
        factor = x*(1-x^2)^(m/2)
        dS = (-1)^m * (factor * dS + dfactor * S)
        S =  (-1)^m * factor * S
    end
    s = iseven(m) ? 1 : -1
    return s*S, s*dS
end
