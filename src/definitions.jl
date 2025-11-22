"""
    oblate_angular_leg(m, n, c, x)
    oblate_angular_leg(m, n, c, λ, x)
    oblate_angular_leg(m, n, c, λ, c2k, x)


"""
function prolate_angular_leg(m, n, c, x)
    λ =  find_eigenvalue(m, n, c)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_ang6_1(m, n, dr, x)
end

function prolate_angular_leg(m, n, c, λ, x)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_ang6_1(m, n, dr, x)
end

function prolate_angular_leg(m, n, c, λ, dr, x)
    return spheroidal_ang6_1(m, n, dr, x)
end

"""
    prolate_angular_leg(m, n, c, x)
    prolate_angular_leg(m, n, c, λ, x)
    prolate_angular_leg(m, n, c, λ, c2k, x)


"""
oblate_angular_leg(m, n, c, x) = prolate_angular_leg(m, n, im*c, x)
oblate_angular_leg(m, n, c, λ, x) = prolate_angular_leg(m, n, im*c, λ, x)
oblate_angular_leg(m, n, c, λ, c2k, x) = prolate_angular_leg(m, n, im*c, λ, c2k, x)


"""
    prolate_angular_ps(m, n, c, x)
    prolate_angular_ps(m, n, c, λ, x)
    prolate_angular_ps(m, n, c, λ, c2k, x)


"""
function prolate_angular_ps(m, n, c, x)
    λ =  find_eigenvalue(m, n, c)
    dr = compute_dr2_mix(m, n, c, λ)
    c2k = compute_c2k(m, n, dr)
    return spheroidal_ang6_2(m, n, c2k, x)
end

function prolate_angular_ps(m, n, c, λ, x)
    dr = compute_dr2_mix(m, n, c, λ)
    c2k = compute_c2k(m, n, dr)
    return spheroidal_ang6_2(m, n, c2k, x)
end

function prolate_angular_ps(m, n, c, λ, c2k, x)
    return spheroidal_ang6_2(m, n, c2k, x)
end

"""
    prolate_angular_ps(m, n, c, x)
    prolate_angular_ps(m, n, c, λ, x)
    prolate_angular_ps(m, n, c, λ, c2k, x)


"""
oblate_angular_ps(m, n, c, x) = prolate_angular_ps(m, n, im*c, x)
oblate_angular_ps(m, n, c, λ, x) = prolate_angular_ps(m, n, im*c, λ, x)
oblate_angular_ps(m, n, c, λ, c2k, x) = prolate_angular_ps(m, n, im*c, λ, c2k, x)


"""
    prolate_radial1(m, n, c, ξ)
    prolate_radial1(m, n, c, λ, ξ)
    prolate_radial1(m, n, c, λ, dr, ξ)


"""
function prolate_radial1(m, n, c, ξ)
    λ =  find_eigenvalue(m, n, c)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_1(m, n, c, dr, ξ)
end

function prolate_radial1(m, n, c, λ, ξ)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_1(m, n, c, dr, ξ)
end

function prolate_radial1(m, n, c, λ, dr, ξ)
    return spheroidal_rad_1(m, n, c, dr, ξ)
end

"""
    prolate_radial2(m, n, c, ξ)
    prolate_radial2(m, n, c, λ, ξ)
    prolate_radial2(m, n, c, λ, dr, ξ)


"""
function prolate_radial2(m, n, c, ξ)
    λ =  find_eigenvalue(m, n, c)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_2(m, n, c, dr, ξ)
end

function prolate_radial2(m, n, c, λ, ξ)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_2(m, n, c, dr, ξ)
end

function prolate_radial2(m, n, c, λ, dr, ξ)
    return spheroidal_rad_2(m, n, c, dr, ξ)
end



"""
    oblate_radial1(m, n, c, ξ)
    oblate_radial1(m, n, c, λ, ξ)
    oblate_radial1(m, n, c, λ, dr, ξ)


"""
function oblate_radial1(m, n, c, ξ)
    λ =  find_eigenvalue(m, n, im*c)
    dr = compute_dr2_mix(m, n, im*c, λ)
    return spheroidal_rad_1(m, n, im*c, dr, im*ξ)
end

function oblate_radial1(m, n, c, λ, ξ)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_1(m, n, im*c, dr, im*ξ)
end

function oblate_radial1(m, n, c, λ, dr, ξ)
    return spheroidal_rad_1(m, n, im*c, dr, im*ξ)
end

"""
    oblate_radial2(m, n, c, ξ)
    oblate_radial2(m, n, c, λ, ξ)
    oblate_radial2(m, n, c, λ, dr, ξ)


"""
function oblate_radial2(m, n, c, ξ)
    λ =  find_eigenvalue(m, n, im*c)
    dr = compute_dr2_mix(m, n, im*c, λ)
    return spheroidal_rad_2(m, n, im*c, dr, im*ξ)
end

function oblate_radial2(m, n, c, λ, ξ)
    dr = compute_dr2_mix(m, n, c, λ)
    return spheroidal_rad_2(m, n, im*c, dr, im*ξ)
end

function oblate_radial2(m, n, c, λ, dr, ξ)
    return spheroidal_rad_2(m, n, im*c, dr, im*ξ)
end


"""
    prolate_cv(m, n, c)


"""
function prolate_cv(m, n, c)
    return find_eigenvalue(m, n, c)
end

"""
    oblate_cv(m, n, c)


"""
function oblate_cv(m, n, c)
    return find_eigenvalue(m, n, im*c)
end

"""
    prolate_cv_seq(m, n, c)


"""
function prolate_cv_seq(m, n, c)
    return find_eigenvalue_seq(m, n, c)
end

"""
    oblate_cv_seq(m, n, c)


"""
function oblate_cv_seq(m, n, c)
    return find_eigenvalue_seq(m, n, im*c)
end
