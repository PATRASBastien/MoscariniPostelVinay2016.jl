# Moscarini, Postel-Vinay (2016) Replication

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PATRASBastien.github.io/MoscariniPostelVinay_Replication.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PATRASBastien.github.io/MoscariniPostelVinay_Replication.jl/dev/)
[![Build Status](https://github.com/PATRASBastien/MoscariniPostelVinay_Replication.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/PATRASBastien/MoscariniPostelVinay_Replication.jl/actions/workflows/CI.yml?query=branch%3Amaster)

This repository contains the Julia package replicating the computational part of the paper:

Moscarini, G., & Postel-Vinay, F. (2016). Wage posting and business cycles: A quantitative exploration. *Review of Economic Dynamics*, 19, 135-160. (https://doi.org/10.1016/j.red.2015.11.001)

This paper presents a new version of the [Burdett and Mortensen (1998)](https://doi.org/10.2307/2527292) model, known as the Stochastic BM (SBM) model, providing a quantitative exploration of business cycles in a frictional labor market under contract-posting. It proposes an algorithm to simulate the equilibrium of this model. 

## Contents

- [**calibration.jl**](calibration.jl): This file contains the `setup_simulation` function, which is used to calibrate the model's parameters to gauge the business cycle properties.
- [**sim_values.jl**](sim_values.jl): This file provides the function `simvalues`, which outputs the structure containing all simulated series.
- [**wagepost_dyn.jl**](wagepost_dyn.jl): This file contains the `simulateModel` function, which ties everything together. It takes the baseline calibration from the paper through `setup_simulation` and runs the simulation using `simvalues`.

## Requirements
- The Julia code was written in v1.8
- **Packages required**:
Random, Distributions, QuantEcon, LinearAlgebra, DelimitedFiles, Statistics, Kronecker, JLD, Plots

## How to Run

We recommend running these lines to compute the figure from the paper:

<pre><code>using MoscariniPostelVinay_Replication
S, P, TB = setup_simulation(291198, 101, 12*5, 12*50, 12*70, .94, .006, 20, 0, 10, 2.5, .95^(1/12), 1, 49, 44, .13)
SIM, S, EXO = simulateModel(0, "noexit")
p = plot_final(SIM, S, P, EXO)
</code></pre>

The output is a figure showing two different versions of the Beveridge curve produced by the model. The blue plot is a plot of the unemployment rate against aggregate job adverts, while the red plot is a plot of aggregate worker search effort against aggregate job adverts.





