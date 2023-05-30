using Random, Distributions, QuantEcon, LinearAlgebra, DelimitedFiles, Statistics, Kronecker


function setup_simulation(seed, gs, trim_beg, trim_end, tmax, rho, sigma, ns, b, pmax, a, beta, m_exp, mcv_exp, mcv_scale, s)

    # SETTINGS (structure S)
    S = Dict()
    # number of aggregate states
    S["ns"] = ns
    # seed
    S["seed"] = seed
    Random.seed!(S["seed"])

    # grid size for p (or Gamma)
    S["gs"] = gs
    # burn-in periods (actual simulation will be for S.tmax - S.trim_beg - S.trim_end periods).
    S["trim_beg"] = trim_beg
    S["trim_end"] = trim_end
    # simulation horizon
    S["tmax"] = tmax + S["trim_beg"] + S["trim_end"]

    # EXOGENOUS PARAMETERS (structure P)
    P = Dict()
    # aggregate productivity shock
    P["rho"] = rho
    P["sigma"] = sigma
    mc = QuantEcon.tauchen(S["ns"], P["rho"], P["sigma"], 0, 4)
    logomega = mc.state_values
    P["trmat"] = mc.p
    P["omega"] = 100 .* exp.(logomega)
    P["sum_trmat"] = cumsum(P["trmat"], dims = 2) - P["trmat"]
    # ergodic distribution of omega
    P["erg_d_om"] = nullspace(transpose(P["trmat"]) .- Matrix(I, size(P["trmat"])))
    P["erg_d_om"] .= P["erg_d_om"] ./ sum(P["erg_d_om"])
    # unemployment income
    P["b"] = b .* ones(size(P["omega"]))
    # job destruction
    EU_rate = DelimitedFiles.readdlm("EU_rate.txt")
    P["delta"] = (maximum(logomega) .- logomega).^2.5
    P["delta"] = (P["delta"] .- (P["erg_d_om"]' * P["delta"])) ./ sqrt(sum((P["erg_d_om"]' .* P["delta"]).^2) .- sum(P["erg_d_om"]' * P["delta"]).^2)
    P["delta"] = mean(EU_rate) .+ std(EU_rate) .* P["delta"]
    # firm type sampling distribution ( ~ Pareto over [1, P.pmax] ):
    P["pmax"] = pmax
    P["a"] = a
    P["ftypepdf"] = x -> P["a"] * x .^ (-P["a"]-1) / (1 - P["pmax"]^(-P["a"]))
    P["ftypecdf"] = x -> (1 - x .^ (-P["a"])) / (1 - P["pmax"]^(-P["a"]))
    P["ftypecdf_inv"] = x -> (1 - (1 - P["pmax"]^(-P["a"])) * x) .^ (-1/P["a"])
    # discount rate
    P["beta"] = beta
    # matching function (Cobb-Douglas)
    P["m_exp"] = m_exp
    P["contact"] = x -> min.(1, x.^P["m_exp"])
    # (inverse of) marginal cost of vacancies
    P["mcv_exp"] = mcv_exp
    P["mcv_scale"] = mcv_scale
    P["mcv_inv"] = x -> x .^ (1/P["mcv_exp"]) / P["mcv_scale"] ^ (1 + 1/P["mcv_exp"])
    # Marginal cost of vacancies
    P["cv"] = x -> (P["mcv_scale"] * x) .^ (1 + P["mcv_exp"]) / (1 + P["mcv_exp"])
    # Worker search effort
    P["s"] = s
    P["seffort"] = (x, k) -> 1 .- x .+ P["s"] .* (1 .- P["delta"][k]) .* x


    ## TOOLBOX (structure TB)
    TB = Dict()
    # nodes for tabulation, integration and interpolation...
    # ... p's
    TB["px"] = P["ftypecdf_inv"].(range(0, stop=.9975, length=(S["gs"] - 9)))
    TB["px"] = [TB["px"][1:end-1]; range(TB["px"][end], stop=P["pmax"], length=10)]

    # ... quantiles
    TB["x"] = P["ftypecdf"].(TB["px"])

    # flow value of a match
    TB["flval"] = kronecker(P["omega"], TB["px"]) .- kronecker(P["b"], ones(S["gs"]))

    # integration operator (trapezoidal rule)
    TB["intop"] = 2 * LowerTriangular(ones(S["gs"], S["gs"])) - Matrix(I, S["gs"], S["gs"])
    TB["intop"][:, 1] .= 1
    TB["intop"][1, 1] = 0
    TB["intop"] = TB["intop"] .* [0; diff(TB["px"])]' ./ 2

    # integration error correction (cheap patch):
    # firm type density integrates to CDF exactly at grid points
    TB["intop"] = TB["intop"] ./ (TB["intop"] * P["ftypepdf"](TB["px"]))
    TB["intop"] = TB["intop"] .* TB["x"]
    TB["intop"][1, :] .= 0

    ## EXPANDED PARAMETERS
    TB["intop_exp"] = kronecker(Matrix(I, S["ns"], S["ns"]), TB["intop"])
    P["delta_exp"] = kronecker(P["delta"], ones(S["gs"]))
    P["trmat_exp"] = kronecker(P["trmat"], I(S["gs"]))
    P["ftypepdf_exp"] = kronecker(ones(S["ns"]), P["ftypepdf"].(TB["px"]))

    return S, P, TB
end