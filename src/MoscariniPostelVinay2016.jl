module MoscariniPostelVinay2016

    using Pkg

    Pkg.instantiate()

    include("calibration.jl")
    include("wagepost_dyn.jl")
    include("sim_values.jl")
    include("plot_final.jl")

end
