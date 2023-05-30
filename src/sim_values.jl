include("calibration.jl")

using Kronecker

@doc """
Returning the structure SIM of all simulated series, including an updated path of job ads a_t(p), 
given an initial guess a0_t(p) and a realization of aggregate states.
...
# Arguments
- `om::int`: omega
- `S::Dict`: Output from calibration.jl
- `TB::Dict`: Output from calibration.jl
- `exit_chk::string`: Check for wether we allow for exit or not
- `a0_in::Vector{Float64}`: sampling distribution for current guess of a(p) 
- `x::Vector{Float64}`: initial guess for x=[0, 0.005, 0.07]
- `P::Dict`: Output from calibration.jl
- `c_st_i::Matrix{Int64}`: aggregate states
- `a0_type::string`: a0_in type
...
"""->
function simvalues(om, c_st_i, S, P, TB, exit_chk, exit, a0_in, x, a0_type)
    # ALLOCATE MEMORY
    # check input for creation of additional variables

    sim = Dict()

    if exit_chk == "exit"
        sim["S"] = zeros(S["ns"] * S["gs"], S["tmax"])
    else
        exit = zeros(S["ns"] * S["gs"], S["tmax"])
    end
    
    if isempty(exit)
        exit = zeros(S["ns"] * S["gs"], S["tmax"])
    end
    
    # allocate rest of memory 

    sim["mu"] = zeros(S["ns"] * S["gs"], S["tmax"])
    sim["pi"] = zeros(S["ns"] * S["gs"], S["tmax"]) 
    sim["V"] = zeros(S["ns"] * S["gs"], S["tmax"])
    sim["a"] = zeros(S["ns"] * S["gs"], S["tmax"])
    sim["dF"] = zeros(S["ns"] * S["gs"], S["tmax"])
    sim["F"] = zeros(S["ns"] * S["gs"], S["tmax"])
    sim["exit"] = zeros(S["ns"] * S["gs"], S["tmax"])
    sim["N"] = zeros(S["ns"] * S["gs"], S["tmax"])
    sim["dN"] = zeros(S["ns"] * S["gs"], S["tmax"])

    sim["lambda"] = zeros(S["ns"], S["tmax"])
    sim["A"] = zeros(S["ns"], S["tmax"])
    sim["intV"] = zeros(S["ns"], S["tmax"])

    lambda_exp = zeros(S["ns"] * S["gs"], S["tmax"]) 
    
    # CONSTRUCT SAMPLING DISTRIBUTION FOR CURRENT GUESS OF a(p)
    if a0_type == "par"
        # guess of a0 given current parameters:
        # approximate a using an affine function of p
        a0_out = x[1] .* kronecker(ones(S["ns"], 1), TB["px"]) .+ x[2]^2 .+ x[3]^2 .* (kronecker(P["omega"], ones(S["gs"], 1)) .- P["omega"][1])
        a0_out = a0_out .* ones(1, S["tmax"])
        a0_out[exit .== 1] .= 0
    else
        a0_out = a0_in
    end 
    
    # implied aggregate vacancies and sampling distribution
    sim["dF"] = a0_out .* P["ftypepdf_exp"]
    sim["A"] = TB["intop_exp"][S["gs"]:S["gs"]:S["ns"]*S["gs"], :] * sim["dF"]
    A_exp = kronecker(sim["A"], ones(S["gs"], 1))
    sim["dF"] = sim["dF"] ./ A_exp
    sim["F"] = TB["intop_exp"] * sim["dF"]

    # SIMULATION OF CONTACT RATES AND SIZE DISTRIBUTION
    # INITIALIZATION (PERIOD 1)
    u_ini = 0.05 
    sim["lambda"][:, 1] = P["contact"](sim["A"][:, 1] ./ P["seffort"](1 - u_ini, 1:S["ns"]))
    lambda_exp[:, 1] = kronecker(sim["lambda"][:, 1], ones(S["gs"], 1))
    sim["dN"][:, 1] = P["delta_exp"] .+ P["s"] .* (1 .- P["delta_exp"]) .* lambda_exp[:, 1] .* (1 .- sim["F"][:, 1])
    sim["dN"][:, 1] = sim["dN"][:, 1] .* sim["dN"][:, 1]
    sim["dN"][:, 1] = (1 - u_ini) .* P["delta_exp"] .* (P["delta_exp"] .+ P["s"] .* lambda_exp[:, 1] .* (1 - P["delta_exp"][om[1]])) ./ sim["dN"][:, 1]
    sim["dN"][:, 1] = sim["dN"][:, 1] .* sim["dF"][:, 1]
    sim["N"][:, 1] = TB["intop_exp"] * sim["dN"][:, 1]


    for t = 2:S["tmax"]
        sim["lambda"][:, t] = P["seffort"](sim["N"][c_st_i[S["gs"], t - 1], t - 1], 1:S["ns"])
        sim["lambda"][:, t] = P["contact"](sim["A"][:, t] ./ sim["lambda"][:, t])
        lambda_exp[:, t] = kronecker(sim["lambda"][:, t], ones(S["gs"], 1))
        sim["dN"][:, t] = (1 .- P["delta_exp"]) .* (1 .- P["s"] * lambda_exp[:, t] .* (1 .- sim["F"][:, t])) .* kronecker(ones(S["ns"], 1), sim["dN"][c_st_i[:, t - 1], t - 1]) + (P["s"] * (1 .- P["delta_exp"]) .* lambda_exp[:, t] .* kronecker(ones(S["ns"], 1), sim["N"][c_st_i[:, t - 1], t - 1]) + lambda_exp[:, t] .* (1 - sim["N"][c_st_i[S["gs"], t - 1], t - 1])) .* sim["dF"][:, t]
        sim["dN"][:, t] = sim["dN"][:, t] .* (1 .- exit[:, t])
        sim["N"][:, t] = TB["intop_exp"] * sim["dN"][:, t]
    end 
    
    # construct hires for substitution in FOCs
    # note: hires at initial date inaccurate, but never used
    H_0 = lambda_exp .* (1 .- kronecker(hcat(sim["N"][S["gs"]:S["gs"]:S["ns"]*S["gs"], 1], sim["N"][S["gs"]:S["gs"]:S["ns"]*S["gs"], 1:(S["tmax"] - 1)]), ones(S["gs"], 1)))
    H_1 = P["s"] .* lambda_exp .* (1 .- P["delta_exp"] * ones(1, S["tmax"]))
    H_2 = hcat(sim["dN"][:, 1], sim["dN"][:, 1:(S["tmax"] - 1)])
    sim["dH"] = H_1 .* H_2
    sim["H"] = TB["intop_exp"] * sim["dH"]
    sim["H"] = H_0 .+ sim["H"]

    
    # SOLVE FOR μ_t and π_t GIVEN PATH OF SAMPLING WEIGHTS
    # FINAL-PERIOD APPROXIMATION: μ_T - U_T = (p - b)
    sim["mu"][:, S["tmax"]] = TB["flval"]
    # if exit is to be checked for, need to compute level of firm value, too
    if exit_chk == "exit"
        sim["S"][:, S["tmax"]] = TB["flval"] .* sim["dN"][:, S["tmax"]]
    end 

    # NOW SOLVE BACKWARD FOR μ_t, V_t, and π_t in earlier periods
    for t = (S["tmax"] - 1):-1:1
        # integral term (construct V(p)-U in passing)
        rhs_int = sim["dH"][:, t + 1] .* sim["mu"][:, t + 1] .* (1 .- exit[:, t + 1])
        rhs_int = TB["intop_exp"] * rhs_int
        sim["V"][:, t + 1] = rhs_int ./ sim["H"][:, t + 1]
        # integral of V
        rhs_int = TB["intop_exp"] * (sim["V"][:, t + 1] .* sim["dF"][:, t + 1])
        sim["intV"][:, t + 1] = rhs_int[S["gs"]:S["gs"]:S["ns"]*S["gs"]]

        # continuation values, state by state
        cval_mu = (1 .- P["delta_exp"]) .* (
            (1 .- P["s"] * lambda_exp[:, t + 1] .* (1 .- sim["F"][:, t + 1])) .* sim["mu"][:, t + 1] .- P["s"] * lambda_exp[:, t + 1] .* rhs_int
        )
        cval_mu = cval_mu .- (1 .- P["s"] * (1 .- P["delta_exp"])) .* lambda_exp[:, t + 1] .* kronecker(sim["intV"][:, t + 1], ones(S["gs"], 1))
        sim["mu"][:, t] = TB["flval"] .+ P["beta"] .* P["trmat_exp"] * cval_mu
        sim["mu"][:, t] = max.(sim["mu"][:, t], 0)
    end  

    # finally construct π_t  
    sim["pi"] = max.(sim["mu"] .- sim["V"], 0)
    sim["exit"] = sim["pi"] .== 0
    # fill in date 1 (it doesn't really matter what with)
    sim["pi"][:, 1] = sim["pi"][:, 2]
    

    ## IMPLIED PATH OF SAMPLING WEIGHTS
    # apply optimality condition
    sim["a"] = P["mcv_inv"](sim["pi"])
    sim["a"] = sim["a"] .* A_exp ./ sim["H"]
    
    dist = sim["a"] .- a0_out
    dist = dist[:, (1 + S["trim_beg"]):(S["tmax"] - S["trim_end"])]
    dist = sqrt(mean(mean(dist .* dist)))
    return dist, sim
end