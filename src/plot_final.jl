
@doc """
plot_final(SIM, S, P, EXO)

Reproduce the plots from the paper to check our results: https://www.aeaweb.org/articles?id=10.1257/aer.p20161051 
...
# Arguments
- `SIM::Dict`: Output from wagepost_gyn.jl
- `S::Dict`: Output from calibration.jl
- `P::Dict`: Output from calibration.jl
- `EXO::Dict`: Output from wagepost_gyn.jl
...
""" ->
function plot_final(SIM, S, P, EXO)
 
    # First plot
    x = 1 .- SIM["N"][S["gs"], :]
    p = plot(x, SIM["A"][1,:], label = "unemployed only")

    # Compute jobseekers
    jobseekers  = 1 .- SIM["N"][S["gs"], :] .+ P["s"] .* (1 .- P["delta"][EXO["omega"]])' .* SIM["N"][S["gs"], :]
    jobseekers  = jobseekers .- mean(jobseekers) .+ mean(1 .- SIM["N"][S["gs"], :])

    # Second plot
    plot!(p, jobseekers, SIM["A"][1,:], label="all job seekers")

    # Labels
    xlabel!("unemployment / total search effort")
    ylabel!("job adverts A_t")

    display(p)

    savefig(p, "figure.png")

    return p
end

