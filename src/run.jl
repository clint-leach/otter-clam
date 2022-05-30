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

# Number of sites
N = size(sag, 3)

# Number of years
T = size(sag, 2)

# Reading in and scaling covariates
X = deepcopy(data[:X])
X_all = deepcopy(data[:X_all])

Xmeans = mean(X_all, dims = 1)
Xsd = std(X_all, dims = 1)

X = (X .- Xmeans) ./ Xsd
X_all = (X_all .- Xmeans) ./ Xsd

# Buidling interpolator to give the solver
λ = [interpolate(λ_grid[:, i], BSpline(Linear())) for i in 1:N]

# Prediction processing
M = size(data[:lambda_all], 2)
λ_all = [interpolate(data[:lambda_all][:, i], BSpline(Linear())) for i in 1:M]

# Setting up model and parameter objects
m = model(z = sag,
          λ = λ,
		  X_r = hcat(fill(1.0, N), X[:, 1:2], X[:, 2] .^ 2),
		  X_a = hcat(fill(1.0, N), X[:, 2]),
		  κ_prior = Gamma(5, 6),
		  ν_prior = Gamma(5.0, 0.02),
		  σ_prior = Beta(4, 1), 
		  Ω_β_r = PDiagMat([0.11, 0.11, 0.11, 0.11]),
		  μ_β_r = [0.0, 0.0, 0.0, 0.0],
		  Ω_β_a = PDiagMat([0.01, 0.01]),
		  μ_β_a = [5.0, 0.0],
		  Doo = data[:Doo],
		  Duu = data[:Duu],
		  Duo = data[:Duo],
		  σ_a = 5.0,
		  ρ_a = 1000.0,
		  σ_r = 3.0,
		  ρ_r = 1000.0,
		  X_r_all = hcat(fill(1.0, M), X_all[:, 1:2], X_all[:, 2] .^ 2),
		  X_a_all = hcat(fill(1.0, M), X_all[:, 2]),
		  λ_all = λ_all,
		  κ_tune = 1.0,
		  ν_tune = 0.002,
		  σ_tune = 0.02,
		  a_tune = ScalMat(N, 0.1),
		  r_tune = ScalMat(N, 5e-4),
		  )

pars = parameters(r = fill(1.0, m.N),
				  β_r = [0.0, 0.0, 0.0, 0.0],
				  a = fill(5.0, m.N),
				  β_a = [0.0, 0.0, 0.0, 0.0],
				  ν = 0.25,
				  κ = 10.0,
				  σ = 0.5, 
				  u = fill(0.0, m.T, m.N),
				  z = fill(0.0, m.T, m.N),
				  loglik = 0.0)

# Running MCMC chain
chain = mcmc(m, pars, 10, 250000, 250000)

R"chain <- readRDS('output/chain.rds')"
@rget chain

# Generating predictions over whole nearshore
preds =  predict(chain, m)

# Saving output to make plots in R
R"saveRDS($(chain), 'output/chain.rds')"
R"saveRDS($(preds), 'output/predictions.rds')"

# Checking acceptance
sum(chain[:accept_kappa])
sum(chain[:accept_nu]) 
sum(chain[:accept_sigma])
sum(chain[:accept_r])
sum(chain[:accept_a])

# Plotting β_r
histogram(chain[:beta_r][1, :], normalize = :true)
plot!(Normal(0.0, 10.0))

histogram(chain[:beta_r][2, :], normalize = :true)
plot!(Normal(0.0, 10.0))

histogram(chain[:beta_r][3, :], normalize = :true)
plot!(Normal(0.0, 10.0))

histogram(chain[:beta_r][4, :], normalize = :true)
plot!(Normal(0.0, 10.0))

# Plotting β_a
histogram(chain[:beta_a][1, :], normalize = :true)
plot!(Normal(5.0, 10.0))

histogram(chain[:beta_a][2, :], normalize = :true)
plot!(Normal(0.0, 10.0))

histogram(chain[:beta_a][3, :], normalize = :true)
plot!(Normal(0.0, 10.0))

# Univariate plots
plot(chain[:nu])
histogram(chain[:nu], normalize = :true)
plot!(m.ν_prior)

plot(chain[:sigma])
histogram(chain[:sigma], normalize = :true)
plot!(m.σ_prior)

plot(chain[:kappa])
histogram(chain[:kappa], normalize = :true)
plot!(m.κ_prior)

# Dynamics at select sites

# Strawberry 
plot(1:26, chain[:u][:, 47, (end-100):end], color = :gray, legend = false)

# Puffin
plot(1:26, chain[:u][:, 42, (end-1000):end], color = :gray, legend = false)

# 64
plot(1:26, chain[:u][:, 22, (end-100):end], color = :gray, legend = false)