using SpheroidalWaveFunctions
using CondaPkg
CondaPkg.add("scipy")
using PythonCall
scipy = pyimport("scipy")

using Test

## CHARACTERISTIC VALUE

@testset "Characteristic value of prolate spheroidal function" begin 
    plus_n = 15
    for c in range(0.01, 101, 31)
        for m = 0:30
            for n in m:(m+plus_n)
                py_cv = pyconvert(Float64, scipy.special.pro_cv(m, n, c))
                @test isapprox(py_cv, prolate_cv(m, n, c))
            end
        end
    end
end

@testset "Characteristic value of oblate spheroidal function" begin 
    plus_n = 15
    for c in range(0.01, 101, 31)
        for m = 0:30
            for n in m:(m+plus_n)
                py_cv = pyconvert(Float64, scipy.special.obl_cv(m, n, c))
                @test isapprox(py_cv, oblate_cv(m, n, c))
            end
        end
    end
end

@testset "Prolate angular function (Legendre)" begin 
    plus_n = 10
    for xi in range(-0.98, 0.98, 10)
        for c in range(0.01, 20.0, 10)
            for m = 0:20
                for n in m:(m+plus_n)
                    py_ang1, py_dang1 = pyconvert(Tuple{Float64, Float64}, scipy.special.pro_ang1(m, n, c, xi))
                    ju_ang1, ju_dang1 = prolate_angular_leg(m, n, c, xi)
                    @test abs((py_ang1 - ju_ang1)/ py_ang1) <= 1e-4 
                    @test abs((py_dang1 - ju_dang1)/ py_dang1) <= 1e-4 
                end
            end
        end
    end
end

@testset "Prolate angular function (Power series)" begin 
    plus_n = 10
    for xi in range(-0.98, 0.98, 10)
        for c in range(0.01, 20.0, 10)
            for m = 0:20
                for n in m:(m+plus_n)
                    py_ang1, py_dang1 = pyconvert(Tuple{Float64, Float64}, scipy.special.pro_ang1(m, n, c, xi))
                    ju_ang1, ju_dang1 = prolate_angular_ps(m, n, c, xi)
                    @test abs((py_ang1 - ju_ang1)/ py_ang1) <= 1e-4 
                    @test abs((py_dang1 - ju_dang1)/ py_dang1) <= 1e-4 
                end
            end
        end
    end
end


@testset "Oblate angular function (Legendre)" begin 
    plus_n = 10
    for xi in range(-0.98, 0.98, 10)
        for c in range(0.01, 15.0, 10)
            for m = 0:20
                for n in m:(m+plus_n)
                    py_ang1, py_dang1 = pyconvert(Tuple{Float64, Float64}, scipy.special.obl_ang1(m, n, c, xi))
                    ju_ang1, ju_dang1 = oblate_angular_leg(m, n, c, xi)
                    @test abs((py_ang1 - ju_ang1)/ py_ang1) <= 1e-4 
                    @test abs((py_dang1 - ju_dang1)/ py_dang1) <= 1e-4 
                end
            end
        end
    end
end

@testset "Oblate angular function (Power series)" begin 
    plus_n = 10
    for xi in range(-0.98, 0.98, 10)
        for c in range(0.01, 15.0, 10)
            for m = 0:20
                for n in m:(m+plus_n)
                    py_ang1, py_dang1 = pyconvert(Tuple{Float64, Float64}, scipy.special.obl_ang1(m, n, c, xi))
                    ju_ang1, ju_dang1 = oblate_angular_ps(m, n, c, xi)
                    @test abs((py_ang1 - ju_ang1)/ py_ang1) <= 1e-4 
                    @test abs((py_dang1 - ju_dang1)/ py_dang1) <= 1e-4 
                end
            end
        end
    end
end

@testset "Prolate Radial function 1" begin 
    plus_n = 10
    for xi in range(1.01, 8.0, 6)
        for c in range(0.01, 30.0, 7)
            for m = 0:1:20
                for n in m:(m+plus_n)
                    py_rad1, py_drad1 = pyconvert(Tuple{Float64, Float64}, scipy.special.pro_rad1(m, n, c, xi))
                    ju_rad1, ju_drad1 = prolate_radial1(m, n, c, xi)
                    @test abs((py_rad1 - ju_rad1)/ py_rad1) <= 1e-4 
                    @test abs((py_drad1 - ju_drad1)/ py_drad1) <= 1e-4 
                end
            end
        end
    end
end

@testset "Prolate Radial function 2" begin 
    plus_n = 10
    for xi in range(1.01, 8.0, 6)
        for c in range(0.01, 30.0, 7)
            for m = 0:1:20
                for n in m:(m+plus_n)
                    py_rad2, py_drad2 = pyconvert(Tuple{Float64, Float64}, scipy.special.pro_rad2(m, n, c, xi))
                    ju_rad2, ju_drad2 = prolate_radial2(m, n, c, xi)
                    @test abs((py_rad2 - ju_rad2)/ py_rad2) <= 1e-4 
                    @test abs((py_drad2 - ju_drad2)/ py_drad2) <= 1e-4 
                end
            end
        end
    end
end


@testset "Oblate Radial function 1" begin 
    plus_n = 10
    for xi in range(1.01, 8.0, 6)
        for c in range(0.01, 30.0, 7)
            for m = 0:1:20
                for n in m:(m+plus_n)
                    py_rad1, py_drad1 = pyconvert(Tuple{Float64, Float64}, scipy.special.obl_rad1(m, n, c, xi))
                    ju_rad1, ju_drad1 = oblate_radial1(m, n, c, xi)
                    @test abs((py_rad1 - ju_rad1)/ py_rad1) <= 1e-4 
                    @test abs((py_drad1 - ju_drad1)/ py_drad1) <= 1e-4 
                end
            end
        end
    end
end

@testset "Oblate Radial function 2" begin 
    plus_n = 10
    for xi in range(1.01, 8.0, 6)
        for c in range(0.01, 30.0, 7)
            for m = 0:1:20
                for n in m:(m+plus_n)
                    py_rad2, py_drad2 = pyconvert(Tuple{Float64, Float64}, scipy.special.obl_rad2(m, n, c, xi))
                    ju_rad2, ju_drad2 = oblate_radial2(m, n, c, xi)
                    @test abs((py_rad2 - ju_rad2)/ py_rad2) <= 1e-4 
                    @test abs((py_drad2 - ju_drad2)/ py_drad2) <= 1e-4 
                end
            end
        end
    end
end