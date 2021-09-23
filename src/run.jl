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
R"data <- readRDS('data/all.rds')"
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
X = deepcopy(data[:X][:, 2:end])
X_all = deepcopy(data[:X_all][:, 2:end])

Xmeans = mean(X, dims = 1)
Xsd = std(X, dims = 1)

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
		  X_r = hcat(fill(1.0, N), X[:, 1:2]),
		  X_a = hcat(fill(1.0, N), X[:, 2], X[:, 2] .^ 2),
		  nq = Int64.(data[:obs_sites][!, :nquad]),
		  κ_prior = Gamma(5, 6),
		  ν_prior = Gamma(1, 0.01),
		  σ_prior = Beta(4, 1), 
		  Ω_β_r = PDiagMat([1.0, 0.25, 0.25]),
		  μ_β_r = [-2.0, 0.0, 0.0],
		  Ω_β_a = PDiagMat([1.0, 0.25, 0.25]),
		  μ_β_a = [3.0, 0.0, 0.0],
		  μ_η_0 = 2.0,
		  Doo = data[:Doo],
		  Duu = data[:Duu],
		  Duo = data[:Duo],
		  σ_a = 0.5,
		  ρ_a = 1000.0,
		  σ_r = 1.0,
		  ρ_r = 1000.0,
		  σ_0 = 1.0,
		  ρ_0 = 1000.0,
		  X_r_all = hcat(fill(1.0, M), X_all[:, 1:2]),
		  X_a_all = hcat(fill(1.0, M), X_all[:, 2], X_all[:, 2] .^ 2),
		  λ_all = λ_all,
		  κ_tune = 12.0,
		  ν_tune = 0.002,
		  σ_tune = 0.02,
		  a_tune = ScalMat(N, 1e-2),
		  r_tune = ScalMat(N, 4e-3),
		  u0_tune = ScalMat(N, 0.04)
		  )

pars = parameters(η_0 = fill(2.0, m.N),
	              r = fill(0.2, m.N),
				  β_r = [0.0, 0.0, 0.0],
				  a = fill(20.0, m.N),
				  β_a = [0.0, 0.0, 0.0],
				  ν = 0.2 / 30.0,
				  κ = 60.0,
				  σ = 0.5, 
				  u = fill(0.0, m.T, m.N),
				  z = fill(0.0, m.T, m.N),
				  loglik = 0.0)

chain = mcmc(m, pars, 1, 0, 1)

preds =  predict(chain, m)

R"saveRDS($(chain), 'output/chain_vara.rds')"
R"saveRDS($(preds), 'output/prediction_nb.rds')"

sum(chain[:accept_kappa])
sum(chain[:accept_nu]) 
sum(chain[:accept_sigma])
sum(chain[:accept_r])
sum(chain[:accept_u0])
sum(chain[:accept_a])

plot(chain[:beta_r][1, :])
plot(chain[:beta_r][2, :])
plot(chain[:beta_r][3, :])

# Univariate plots
plot(chain[:nu])
histogram(chain[:nu], normalize = :true)
plot!(m.ν_prior)

plot(chain[:sigma])
histogram(chain["sigma"], normalize = :true)
plot!(m.σ_prior)

plot(chain[:kappa])
histogram(chain[:kappa], normalize = :true)
plot!(m.κ_prior)

plot(chain[:a][1, :])
histogram(chain[:a], normalize = :true)
plot!(m.a_prior)

# Dynamics

# Strawberry
plot(1:26, chain[:u][:, 47, (end-100):end], color = :gray, legend = false)

# Puffin
plot(1:26, chain[:u][:, 42, (end-100):end], color = :gray, legend = false)

