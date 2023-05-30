module MoscariniPostelVinay2016

    using Pkg, Plots, Random, Distributions, QuantEcon, LinearAlgebra, DelimitedFiles, Statistics, Kronecker, JLD

    Pkg.instantiate()

    include("calibration.jl")
    include("wagepost_dyn.jl")
    include("sim_values.jl")
    include("plot_final.jl")

end
