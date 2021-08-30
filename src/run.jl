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
λ_grid = data[:lambda][6:26, :]

# Reading in prey data
sag = data[:y][:, 6:26, :]

# Number of sites
N = size(sag, 3)

# Number of years
T = size(sag, 2)

# Reading in and scaling covariates
X = deepcopy(data[:X])
X_all = deepcopy(data[:X_all])

Xmeans = mean(X, dims = 1)
Xsd = std(X, dims = 1)

X[:, 2] = (X[:, 2] .- Xmeans[2]) ./ (2 * Xsd[2])
X[:, 3] = (X[:, 3] .- Xmeans[3]) ./ (2 * Xsd[3])

X_all[:, 2] = (X_all[:, 2] .- Xmeans[2]) ./ (2 * Xsd[2])
X_all[:, 3] = (X_all[:, 3] .- Xmeans[3]) ./ (2 * Xsd[3])

# Buidling interpolator to give the solver
λ = [interpolate(λ_grid[:, i], BSpline(Linear())) for i in 1:N]

# Prediction processing
M = size(data[:lambda_all], 2)
λ_all = [interpolate(data[:lambda_all][6:26, i], BSpline(Linear())) for i in 1:M]

# Setting up model and parameter objects
m = model(z = sag,
          λ = λ,
		  X = X,
		  nq = Int64.(data[:obs_sites][!, :nquad]),
		  a_prior = Gamma(2, 10),
		  κ_prior = Gamma(5, 6),
		  K_prior = Gamma(10, 6),
		  σ_prior = Beta(4, 1), 
		  Ω_β_r = PDiagMat([1.0, 0.25, 0.25]),
		  μ_β_r = [-2.0, 0.0, 0.0],
		  Ω_β_0 = PDiagMat([1.0, 0.16, 0.16]),
		  μ_β_0 = [0.0, 0.0, 0.0],
		  Doo = data[:Doo],
		  Duu = data[:Duu],
		  Duo = data[:Duo],
		  σ_r = 1.0,
		  ρ_r = 1000.0,
		  σ_0 = 2.5,
		  ρ_0 = 1000.0,
		  X_all = X_all,
		  λ_all = λ_all,
		  a_tune = 2.0,
		  κ_tune = 3.0,
		  K_tune = 3.0,
		  σ_tune = 0.02,
		  r_tune = ScalMat(N, 5e-3),
		  u0_tune = ScalMat(N, 1e-2)
		  )

pars = parameters(η_0 = fill(0.0, m.N),
				  r = fill(0.2, m.N),
				  β_r = [0.0, 0.0, 0.0],
				  β_0 = [0.0, 0.0, 0.0],
				  a = 20.0,
				  K = 60.0,
				  κ = 60.0,
				  σ = 0.5, 
				  u = fill(0.0, m.T, m.N),
				  z = fill(0.0, m.T, m.N),
				  loglik = fill(0.0, m.N))

chain = mcmc(m, pars, 10, 250000, 250000)

preds =  predict(chain, m)

R"saveRDS($(chain), 'output/chain_nb.rds')"
R"saveRDS($(preds), 'output/prediction_nb.rds')"

sum(chain[:accept_a])
sum(chain[:accept_kappa])
sum(chain[:accept_K]) 
sum(chain[:accept_sigma])
sum(chain[:accept_r]) 
sum(chain[:accept_u0])

plot(chain[:beta_0][1, :])
plot(chain[:beta_0][2, :])
plot(chain[:beta_0][3, :])

plot(chain[:beta_r][1, :])
plot(chain[:beta_r][2, :])
plot(chain[:beta_r][3, :])

# Univariate plots

plot(chain[:K])
histogram(chain[:K], normalize = :true)
plot!(m.K_prior)

plot(chain[:sigma])
histogram(chain["sigma"], normalize = :true)
plot!(m.σ_prior)

plot(chain[:kappa])
histogram(chain[:kappa], normalize = :true)
plot!(m.κ_prior)

plot(chain[:a])
histogram(chain[:a], normalize = :true)
plot!(m.a_prior)

# Dynamics
zmean = mean(chain[:zpred], dims = 3)
plot(1:21, chain[:zpred][:, 39, (end-100):end], color = :gray, legend = false)
plot!(1:21, zmean[:, 39, 1], color = :black, width = 2)
scatter!(1:21, m.z[:, :, 39]', color = :black)