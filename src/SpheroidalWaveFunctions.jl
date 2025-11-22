module SpheroidalWaveFunctions

    using LinearAlgebra: Tridiagonal, SymTridiagonal, eigen
    using Bessels: sphericalbesselj, sphericalbessely

    include("./common.jl")
    include("./characteristic_value.jl")
    include("./coefficients.jl")
    include("./angular_functions.jl")
    include("./radial_functions.jl")
    include("./definitions.jl")

    export prolate_angular_leg, prolate_angular_ps
    export oblate_angular_leg, oblate_angular_ps
    export prolate_radial1
    export prolate_radial2
    export oblate_radial1
    export oblate_radial2

    export prolate_cv
    export oblate_cv

    export prolate_cv_seq
    export oblate_cv_seq

end

