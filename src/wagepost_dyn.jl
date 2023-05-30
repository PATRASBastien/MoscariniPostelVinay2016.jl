# Include the other files
include("calibration.jl")
include("sim_values.jl")

using Random, LinearAlgebra, DelimitedFiles, JLD, Kronecker

function simulateModel(skip_1st_step, exit_chk)
    
    # skip_1st_step: set to 1 to skip first (parametric) step in calibration of vacancies, and to anything else to not skip it
    # exit_chk: "exit" or "noexit", allowing for firm entry and exit

    # set calibrated parameter values
    S, P, TB = setup_simulation(291198, 101, 12*5, 12*50, 12*70, .94, .006, 20, 0, 10, 2.5, .95^(1/12), 1, 49, 44, .13)

    # allocate memory for exogenous draws (structure EXO)
    EXO = Dict("omega" => zeros(Int, 1, S["tmax"]), "c_st_i" => zeros(Int, S["gs"], S["tmax"]))

    ## DRAW SERIES OF AGGREGATE STATES
    dice = rand(1, S["tmax"])
    EXO["omega"][1] = ceil(Int, S["ns"] / 2)

    for t = 2:S["tmax"]
        EXO["omega"][t] = sum(dice[t] .> P["sum_trmat"][EXO["omega"][t-1], :])
        current_state_index = ((EXO["omega"][t] - 1) * S["gs"] + 1):(EXO["omega"][t] * S["gs"])
        EXO["c_st_i"][:, t] = current_state_index[:]
    end

    # fill in current state index at date 1
    current_state_index = ((EXO["omega"][1] - 1) * S["gs"] + 1):(EXO["omega"][1] * S["gs"])
    EXO["c_st_i"][:, 1] = current_state_index[:]

    # Housekeeping
    for varname in names(Main)
        if endswith(string(varname), "_ini") || varname ∈ [:t, :dice, :s, :logomega]
            eval(Meta.parse("$(varname) = nothing"))
        end
    end

    ## SOLVE FOR PATH OF VACANCIES
    # initial guess: optional first parametric round
    if skip_1st_step != 1
        x_a = [0, 0.005, 0.07]
        exit0 = zeros(S["ns"] * S["gs"], S["tmax"])
        _, SIM = simvalues(EXO["omega"], EXO["c_st_i"], S, P, TB, exit_chk, exit0, [], x_a, "par")
        a0 = SIM["a"]
        dist = 10
        save("a_initial_guess.jld", "a0", a0, "exit0", exit0, "dist", dist)
    end

    # second non-parametric round, starting from some initial guess (e.g. a parametric approximation obtained before)
    load("a_initial_guess.jld", "a0", "exit0", "dist")
    i = 0
    relmax = 0.5
    relaxation = relmax

    if dist <= 1e-3
        dist = 10
    end
    println("")

    while dist > 1e-3 && i < 100
        dist_1 = dist
        i += 1
        dist, SIM_cand = simvalues(EXO["omega"], EXO["c_st_i"], S, P, TB, exit_chk, exit0, a0, [], "nonpar")
        println("iteration ", i, ": distance = ", dist, ", relaxation = ", relaxation)
        if dist < dist_1 || i == 1
            a0_1 = a0
            exit0_1 = exit0
            a0 = a0 + relaxation * (SIM_cand["a"] - a0)
            exit0 = SIM_cand["exit"]
            SIM = SIM_cand
        else
            relaxation = 0.5 * relaxation
            if relaxation > relmax / 128
                println("Distance increasing. Reducing relaxation to $relaxation and trying again.")
            else
                println("Distance increasing. Resetting relaxation to $relmax and trying again.")
                relaxation = relmax
            end
            a0 = a0_1 + relaxation * (a0 - a0_1)
            exit0 = exit0_1
            dist = dist_1
        end
    end
    save("a_initial_guess.jld", "a0", a0, "exit0", exit0, "dist", dist)

    # Housekeeping
    for varname in names(Main)
        if startswith(string(varname), "dist") || startswith(string(varname), "a0") || startswith(string(varname), "exit") || varname ∈ [:i, :relaxation, :skip_1st_step, :relmax, :SIM_cand, :x_a]
            eval(Meta.parse("$(varname) = nothing"))
        end
    end

    ## CONSTRUCT OTHER SIMULATED SERIES OF INTEREST
    SIM["U"] = zeros(S["ns"], S["tmax"])
    SIM["w"] = zeros(S["gs"], S["tmax"])

    # time saver
    lambda_exp = kronecker(SIM["lambda"], ones(S["gs"], 1))

    # END-PERIOD APPROXIMATION FOR U_T
    SIM["U"] = P["beta"] * ones(S["ns"], 1) * (P["erg_d_om"]' * (SIM["lambda"] .* SIM["intV"]))
    SIM["U"][:, S["tmax"]] = (P["b"] + mean(SIM["U"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])], dims=2)) / (1 - P["beta"])

    # SOLVE BACKWARD FOR wages and U_t, and check that \mu_t is positive
    for t in (S["tmax"] - 1):-1:1
        cval = SIM["U"][:, t + 1] + SIM["lambda"][:, t + 1] .* SIM["intV"][:, t + 1]
        SIM["U"][:, t] = P["b"] + P["beta"] * P["trmat"] * cval

        cval = (1 .- P["delta_exp"]) .* (1 .- P["s"] * lambda_exp[:, t + 1] .* (1 .- SIM["F"][:, t + 1])) .* SIM["V"][:, t + 1]
        cval .-= (1 .- P["s"] * (1 .- P["delta_exp"])) .* lambda_exp[:, t + 1] .* kronecker(SIM["intV"][:, t + 1], ones(S["gs"], 1))
        cval .-= P["s"] * lambda_exp[:, t + 1] .* (1 .- P["delta_exp"]) .* (TB["intop_exp"] * (SIM["dH"][:, t + 1] .* SIM["mu"][:, t + 1] .* (1 .- SIM["exit"][:, t + 1])) ./ SIM["H"][:, t + 1])
        cval = P["trmat_exp"] * cval
        SIM["w"][:, t] = P["beta"] * cval[EXO["c_st_i"][:, t]]
        SIM["w"][:, t + 1] = SIM["V"][EXO["c_st_i"][:, t + 1], t + 1] .+ P["b"][EXO["omega"][t + 1]] .- SIM["w"][:, t + 1]
    end

    # Housekeeping
    for varname in names(Main)
        if endswith(string(varname), "_exp") || startswith(string(varname), "U_") || varname ∈ [:t, :s, :cval]
            eval(Meta.parse("$(varname) = nothing"))
        end
    end

    delete!(SIM, "dH")
    delete!(SIM, "H")
    delete!(SIM, "intV")

    # Trim to get rid of end-period approx error on U and first-period approx error on N
    EXO["omega"] = EXO["omega"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["N"] = SIM["N"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["dN"] = SIM["dN"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["mu"] = SIM["mu"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["pi"] = SIM["pi"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["V"] = SIM["V"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["U"] = SIM["U"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["w"] = SIM["w"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["a"] = SIM["a"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["A"] = SIM["A"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["F"] = SIM["F"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["dF"] = SIM["dF"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["lambda"] = SIM["lambda"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    SIM["exit"] = SIM["exit"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    EXO["c_st_i"] = EXO["c_st_i"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]

    if haskey(SIM, "S")
        SIM["S"] = SIM["S"][:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    end

    S["tmax"] = S["tmax"] - S["trim_beg"] - S["trim_end"]
    SIM["N_assess"] = SIM["N"]

    return SIM, S, EXO
end