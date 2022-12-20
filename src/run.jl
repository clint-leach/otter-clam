using DifferentialEquations
using Interpolations
using Parameters
using Distributions
using Plots
using StatsPlots
using RCall
using DataFrames
using PDMats
using LinearAlgebra
using StatsFuns

include("model.jl")
include("process.jl")
include("sample.jl")
include("predict.jl")

# Reading in data
R"data <- readRDS('output/all.rds')"
@rget data

# Reading in sea otter forcing time series
λ_grid = data[:lambda]

# Reading in prey data
sag = data[:y]

# Computing number of quadrats sampled
quads = sum(ismissing.(sag) .== 0, dims = 1)[1, :, :]

# Number of sites
N = size(sag, 3)

# Number of years
T = size(sag, 2)

# Reading in and scaling covariates
X = deepcopy(data[:X])

Xmeans = mean(X, dims = 1)
Xsd = std(X, dims = 1)

X_scaled = (X .- Xmeans) ./ Xsd

# Buidling interpolator to give the solver
λ = [interpolate(λ_grid[:, i], BSpline(Linear())) for i in 1:N]

# Building model object with model constants
m = model(z = sag,
          λ = λ,
		  X_r = hcat(fill(1.0, N), X_scaled[:, 1], X_scaled[:, 2], X_scaled[:, 2] .^ 2, X_scaled[:, 3]),
		  X_a = hcat(fill(1.0, N), X_scaled[:, 2], X_scaled[:, 3]),
		  quads = quads,
		  ν_prior = Gamma(5.0, 0.02),
		  σ_prior = Beta(2, 1), 
		  Ω_β_r = PDiagMat(fill(0.25, 5)),
		  μ_β_r = vcat(2.3, fill(0.25, 4)),
		  Ω_β_a = PDiagMat(fill(1.0, 3)),
		  μ_β_a = vcat(-2.0, fill(0.0, 2)),
		  Doo = data[:Doo],
		  σ_a = 1.0,
		  ρ_a = 100.0,
		  σ_r = 1.0,
		  ρ_r = 100.0,
		  ν_tune = 0.01,
		  σ_tune = 0.01,
		  a_tune = ScalMat(N, 0.02),
		  r_tune = ScalMat(N, 0.003)
		  )

# Initializing parameters object
pars = parameters(r = fill(1.0, m.N),
				  β_r = m.μ_β_r,
				  β_end = m.μ_β_r,
				  a = fill(0.1, m.N),
				  β_a = m.μ_β_a,
				  ν = 0.05,
				  σ = fill(0.5, m.N), 
				  u = fill(0.0, m.T, m.N),
				  z = fill(0.0, m.T, m.N),
				  loglik = 0.0)

# Running MCMC chain
chain = mcmc(m, pars, 20, 400000, 400000)

# Saving output to make plots in R
R"saveRDS($(chain), 'output/chain_main.rds')"

# Checking acceptance
sum(chain[:accept_nu]) 
sum(chain[:accept_sigma])
sum(chain[:accept_r])
sum(chain[:accept_a])

# Exploratory Plotting

# Plotting β_r
plot(chain[:beta_r][1, :])
plot!(chain[:beta_end][1, :])

# Plotting β_a
plot(chain[:beta_a][1, :])

# Plotting ν
plot(chain[:nu])

# Plotting σ's
plot(chain[:sigma]', legend = false)

# Dynamics at select sites

# Boulder
plot(1:20, chain[:u][:, 28, (end-100):end], color = :gray, legend = false)