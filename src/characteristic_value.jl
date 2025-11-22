
"""
    build_tridiagonal(m, n, c, N) -> Tridiagonal

Construct the tridiagonal matrix for the three-term recurrence relation of spheroidal wave functions.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter
- `N::Integer`: matrix dimension

# Returns
- Tridiagonal matrix with elements αᵣ, βᵣ, γᵣ
"""
function build_tridiagonal(m, n, c, N)
    start = iseven(m - n) ? 0 : 1
    dd = zeros(N)
    du = zeros(N - 1)
    dl = zeros(N - 1)
    γ² = real(c^2)
    for (idx, r) in enumerate(start:2:start+2*(N-1))
        v1 = (2 * m + r)
        v2 = (2 * m + 2 * r)
        v3 = m + r
        αᵣ = γ² * (v1 + 2) * (v1 + 1) / (v2 + 5) / (v2 + 3)
        βᵣ = v3 * (v3 + 1) + γ² * (2 * v3 * (v3 + 1) - 2 * m^2 - 1) / ((v2 - 1) * (v2 + 3))
        γᵣ = γ² * (r * (r - 1)) / ((v2 - 3) * (v2 - 1))
        dd[idx] = βᵣ
        if idx > 1
            dl[idx-1] = γᵣ
        end
        if idx < N
            du[idx] = αᵣ
        end
    end
    return Tridiagonal(dl, dd, du)
end

"""
    symmetrize_tridiagonal(A) -> SymTridiagonal

Convert a tridiagonal matrix to a symmetric tridiagonal form for eigenvalue computation.

# Arguments
- `A::Tridiagonal`: input tridiagonal matrix

# Returns
- Symmetric tridiagonal matrix with off-diagonal elements √(Aᵢⱼ * Aⱼᵢ)
"""
function symmetrize_tridiagonal(A)    
    n = length(A.d)                                                                                                                                    
    ev  = zeros(Float64, n-1)      # new symmetric diagonal    
    #ev = A.du                    
    for i in 1:n-1                                                          
        ev[i] = sqrt(A.du[i] * A.dl[i])                                              
    end
    return SymTridiagonal(A.d, ev)
end

"""
    cv_matrix(m, n, c, N=10 + div(n - m, 2) + round(Int, abs(c))) -> λ

Compute the characteristic value using the matrix eigenvalue method.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter
- `N::Integer`: matrix dimension (default: adaptive based on n, m, c)

# Returns
- Characteristic value λₘₙ(c)
"""
function cv_matrix(m, n, c, N = 10 + div(n - m, 2) + round(Int, abs(c)))
    p = div(n - m, 2) + 1
    A = build_tridiagonal(m, n, c, N)
    A2 = symmetrize_tridiagonal(A)
    eigv, _ = eigen(A2, 1:p)
    λₘₙ = eigv[p]
    return λₘₙ
end

"""
    cv_matrix_seq(m, n, c, N=10 + div(n - m, 2) + round(Int, abs(c))) -> Vector{λ}

Compute a sequence of characteristic values using the matrix eigenvalue method.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: maximum mode number
- `c::Number`: spheroidal parameter
- `N::Integer`: matrix dimension (default: adaptive based on n, m, c)

# Returns
- Vector of characteristic values
"""
function cv_matrix_seq(m, n, c, N = 10 + div(n - m, 2) + round(Int, abs(c)))
    p = div(n - m, 2) + 1
    A = build_tridiagonal(m, n, c, N)
    A2 = symmetrize_tridiagonal(A)
    eigv, _ = eigen(A2, 1:p)
    λₘₙ = eigv[1:p]
    return λₘₙ
end

"""
    Nᵣᵐ1(m, n, c, λ, Nmax) -> Number

Compute the backward continued fraction for the characteristic value equation.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter
- `λ::Number`: trial characteristic value
- `Nmax::Integer`: maximum number of terms in continued fraction

# Returns
- Value of the backward continued fraction
"""
function Nᵣᵐ1(m, n, c, λ, Nmax)
    γ² = real(c^2)
    v = iseven(n - m) ? 0 : 1
    result = 0.0
    rs = max(n - m - 2 * Nmax, 2 + v):2:(n-m)
    for ri in rs
        result = βᵐ(m, ri, γ²) / (γᵐ(m, ri - 2, γ²) - λ - result)
    end
    return γᵐ(m, n - m, γ²) - λ - result
end

"""
    Nᵣᵐ2(m, n, c, λ, Nmax) -> Number

Compute the forward continued fraction for the characteristic value equation.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter
- `λ::Number`: trial characteristic value
- `Nmax::Integer`: maximum number of terms in continued fraction

# Returns
- Value of the forward continued fraction
"""
function Nᵣᵐ2(m, n, c, λ, Nmax)
    γ² = real(c^2)
    result = 0.0
    rs = (n-m+2*Nmax):-2:(n-m+2)
    for ri in rs
        result = βᵐ(m, ri, γ²) / (γᵐ(m, ri, γ²) - λ - result)
    end
    return result
end


"""
    Nᵣᵐ1_prime(m, n, c, λ, Nmax) -> Number

Compute the derivative of the backward continued fraction with respect to λ.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter
- `λ::Number`: trial characteristic value
- `Nmax::Integer`: maximum number of terms

# Returns
- Derivative dNᵣᵐ1/dλ
"""
function Nᵣᵐ1_prime(m, n, c, λ, Nmax)
    γ² = real(c^2)
    R_prime_in = 0.0 
    R_in = 0.0       
    v = iseven(n - m) ? 0 : 1
    
    rs = max(n - m - 2 * Nmax, 2 + v):2:(n-m)
    
    for ri in rs
        R_next = R_in
        R_prime_next = R_prime_in
        βr = βᵐ(m, ri, γ²) 
        den = γᵐ(m, ri - 2, γ²) - λ - R_next 
        R_in = βr / den
        R_prime_in = (R_in^2 / βr) * (1 + R_prime_next)
    end
    
    return -1.0 - R_prime_in
end

"""
    Nᵣᵐ2_prime(m, n, c, λ, Nmax) -> Number

Compute the derivative of the forward continued fraction with respect to λ.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter
- `λ::Number`: trial characteristic value
- `Nmax::Integer`: maximum number of terms

# Returns
- Derivative dNᵣᵐ2/dλ
"""
function Nᵣᵐ2_prime(m, n, c, λ, Nmax)
    γ² = real(c^2)
    R_prime = 0.0 
    R = 0.0       
    
    rs = (n-m+2*Nmax):-2:(n-m+2)
    for ri in rs
        R_next = R
        R_prime_next = R_prime
        βr = βᵐ(m, ri, γ²)
        den = γᵐ(m, ri, γ²) - λ - R_next 
        R = βr / den
        R_prime = (R^2 / βr) * (1 + R_prime_next)
    end
    
    return R_prime
end

"""
    find_eigenvalue(m, n, c; tol=1e-12, max_iter=25, num_terms=100) -> λ

Compute the characteristic value using Newton-Raphson iteration with continued fractions.

Uses the matrix method to obtain an initial guess, then refines it using Newton-Raphson
iteration on the continued fraction equation.

# Arguments
- `m::Integer`: azimuthal order (m ≥ 0)
- `n::Integer`: mode number (n ≥ m)
- `c::Number`: spheroidal parameter
- `tol::Float64`: convergence tolerance (default: 1e-12)
- `max_iter::Int`: maximum number of iterations (default: 25)
- `num_terms::Int`: number of terms in continued fraction (default: 100)

# Returns
- Characteristic value λₘₙ(c)

# Examples
```julia
λ = find_eigenvalue(1, 2, 0.5)
λ = find_eigenvalue(1, 2, 0.5, tol=1e-14, max_iter=50)
```
"""
function find_eigenvalue(m, n, c; tol = 1e-12, max_iter=25, num_terms=100)
    γ² = real(c^2)
    λ_guess = cv_matrix(m, n, c)
    λ_current = find_eigenvalue_guess(m, n, c, λ_guess; tol = tol, max_iter = max_iter, num_terms = num_terms)
    return λ_current
end

"""
    find_eigenvalue_guess(m, n, c, λ_guess; tol=1e-12, max_iter=25, num_terms=40) -> λ

Refine a characteristic value using Newton-Raphson iteration with a given initial guess.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter
- `λ_guess::Number`: initial guess for the characteristic value
- `tol::Float64`: convergence tolerance (default: 1e-12)
- `max_iter::Int`: maximum number of iterations (default: 25)
- `num_terms::Int`: number of terms in continued fraction (default: 40)

# Returns
- Refined characteristic value λₘₙ(c)
"""
function find_eigenvalue_guess(m, n, c, λ_guess; tol = 1e-12, max_iter=25, num_terms=40)
    λ_current = λ_guess
    for _ in 1:max_iter
        U1 = Nᵣᵐ1(m, n, c, λ_current, num_terms)
        U2 = Nᵣᵐ2(m, n, c, λ_current, num_terms)
        U_val = U1 - U2 
        abs(U_val) <= tol && break
        dU1dλ = Nᵣᵐ1_prime(m, n, c, λ_current, num_terms)
        dU2dλ = Nᵣᵐ2_prime(m, n, c, λ_current, num_terms)
        dUdλ = dU1dλ - dU2dλ
        if abs(dUdλ) < eps(λ_current) 
             break
        end
        λ_new = λ_current - U_val / dUdλ
        λ_current = λ_new
    end
    return λ_current
end

"""
    find_eigenvalue_seq(m, n, c; max_iter=20, num_terms=30) -> Vector{λ}

Compute a sequence of characteristic values for modes m through n.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: maximum mode number
- `c::Number`: spheroidal parameter
- `max_iter::Int`: maximum iterations per eigenvalue (default: 20)
- `num_terms::Int`: number of terms in continued fraction (default: 30)

# Returns
- Vector of characteristic values [λₘₘ, λₘ₊₁,ₘ, ..., λₙₘ]
"""
function find_eigenvalue_seq(m, n, c; max_iter=20, num_terms=30)
    γ² = real(c^2)
    λn1 = cv_matrix_seq(m, n, c)
    λn2 = cv_matrix_seq(m, n - 1, c)
    λn = sort!(vcat(λn1, λn2))
    λ_seq = [find_eigenvalue_guess(m, n, c, λi; max_iter=max_iter, num_terms=num_terms) for (n, λi) in zip((m:m+n), λn)]
    return λ_seq
end


"""
    ps_λ(m, n, c) -> λ

Compute the characteristic value using power series expansion for small c.

Uses the asymptotic expansion λₘₙ(c) ≈ n(n+1) + ℓ₂c² + ℓ₄c⁴ + ℓ₆c⁶ for small values
of the spheroidal parameter c.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter (should be small)

# Returns
- Approximate characteristic value λₘₙ(c)

# Reference
Page 754, Abramowitz and Stegun, "Handbook of Mathematical Functions"
"""
function ps_λ(m, n, c)
    ℓ₀ = n * (n + 1)
    ℓ₂ = 1 / 2 * (-1 - ((2 * m - 1) * (2 * m + 1)) / ((2 * n - 1) * (2 * n + 3)))

    num_4a = (n - m - 1) * (n - m) * (n + m - 1) * (n + m)
    den_4a = (2 * n - 3) * (2 * n - 1)^3 * (2 * n + 1)
    num_4b = (n - m + 1) * (n - m + 2) * (n + m + 1) * (n + m + 2)
    den_4b = (2 * n + 1) * (2 * n + 3)^3 * (2 * n + 5)
    ℓ₄ = 1 / 2 * (num_4a / den_4a - num_4b / den_4b)

    num_6a = num_4b
    den_6a = (2 * n - 1) * (2 * n + 1) * (2 * n + 3)^5 * (2 * n + 5) * (2 * n + 7)

    num_6b = (n - m - 1) * (n - m) * (n + m - 1) * (n + m)
    den_6b = (2 * n - 5) * (2 * n - 3) * (2 * n - 1)^5 * (2 * n + 1) * (2 * n + 3)

    ℓ₆ = (4 * m^2 - 1) * (num_6a / den_6a - num_6b / den_6b)

    v = ℓ₀ * real(c^(2 * 0)) + ℓ₂ * real(c^(2 * 1)) + ℓ₄ * real(c^(2 * 2)) + ℓ₆ * real(c^(2 * 3))
    return v
end


"""
    asymp_λ_pro(m, n, c) -> λ

Compute the characteristic value using asymptotic expansion for large c (oblate case).

Uses the asymptotic expansion for large values of the spheroidal parameter c in the
prolate case.

# Arguments
- `m::Integer`: azimuthal order
- `n::Integer`: mode number
- `c::Number`: spheroidal parameter (should be large)

# Returns
- Approximate characteristic value λₘₙ(c) for large c

# Reference
Page 754, Abramowitz and Stegun, "Handbook of Mathematical Functions"
"""
function asymp_λ_pro(m, n, c)                                                                                         
    q = 2*(n-m)+1                                                                                                                                                                              
    v = c * q + m^2 - 1/8 * (q^2 + 5) - q/(64*c) * (q^2 + 11 - 32^2) - 1/(1024*c^2) * (5*(q^4 + 26*q^2+21) - 384*m^2 * (q^2 + 1)) 
    #v -= 1/c^3 * ( (1/128^2*(33*q^5 + 1594*q^3 + 5621*q) - m^2/128 * (37*q^3 + 167*q) + m^4/8 * q)) -1/c^5 * (1/1024^2 * (527*q^7 + 61529*q^5 + 1043961*q^3 + 2241599*q) - m^2 / (32*1024) * (5379*q^5 + 127550*q^3 + 298951*q) + m^4 / 512 * (355*q^3 + 1505*q) - m^6*q)
    return v                                                                                                      
end