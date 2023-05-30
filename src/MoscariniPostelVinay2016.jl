module MoscariniPostelVinay2016

    using Plots, Random, Distributions, QuantEcon, LinearAlgebra, DelimitedFiles, Statistics, Kronecker

    include("calibration.jl")
    include("wagepost_dyn.jl")
    include("sim_values.jl")
    include("plot_final.jl")

end
