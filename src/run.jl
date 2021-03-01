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

include("model.jl")
include("process.jl")
include("sample.jl")

# Reading in sea otter forcing time series
R"data <- readRDS('data/otter_array.rds')"
@rget data

# Reading in prey data
R"sag <- readRDS('data/SAG_array.rds')"
@rget sag

# Reading in distance matrix
R"dists <- readRDS('data/distances.rds')"
@rget dists

# Subsetting to sites with the prey observed
obs = sum(1 .- ismissing.(sag), dims = 1)[1, :]
obs = findall(obs .> 1)

# obs = [8, 36]
data = data[:, obs]
sag = sag[:, obs]
D = dists[obs, obs]

# Number of sites
N = size(data)[2]

# Number of years
T = size(data)[1]

# Buidling interpolator to give the solver
λ = [interpolate(data[:, i], BSpline(Linear())) for i in 1:N]

# Building fixed covariance matrix for log_r
Σ_r  = PDMat(2.0 * exp.(-0.5 * (D ./ 10000)))
Σ_K  = PDMat(1.0 * exp.(-0.5 * (D ./ 10000)))

# Simulating data
# log_r = rand(MvNormal(fill(log(0.2), N), Σ_r))
# log_K = rand(MvNormal(fill(log(500), N), Σ_K))
#
# p = [log_r, log_K, 50.0, 500.0, λ]
# prob = ODEProblem(prey_all!, rand(Normal(1000, 200), N), (1.0, 26.0), p)
# sol = solve(prob, Tsit5(), saveat=1.0)
#
# z = [rand(truncated(Normal(sol[i, t], 50), 0.0, Inf)) for t in 1:T, i in 1:N]

# Setting up model and parameter objects
m = model(z = sag,
          T = T,
		  tobs = 1:T,
          tspan = (1.0, 26.0),
		  N = N,
          λ = λ,
		  α_r_tune = ScalMat(N, 0.1),
		  α_r_prior = MvNormal(N, 1.0),
		  L_r = cholesky(Σ_r).U,
		  μ_r = fill(log(0.2), N),
		  a_tune = 50.0,
		  a_prior = Gamma(1, 50),
		  κ_tune = 100.0,
		  κ_prior = Gamma(50, 10),
		  α_K_tune = ScalMat(N, 0.1),
		  α_K_prior = MvNormal(N, 1.0),
		  μ_K = fill(log(500.0), N),
		  L_K = cholesky(Σ_K).U,
		  u0_tune = 500.0,
		  u0_prior = Gamma(5, 200),
		  σ_tune = 50,
		  σ_prior = truncated(Normal(0, 100), 0.0, Inf)
		  )

pars = parameters(u0 = fill(1000.0, N),
                  log_r = fill(log(0.2), N),
				  α_r = fill(0.0, N),
				  log_K = fill(log(500.0), N),
				  α_K = fill(0.0, N),
				  a = 50.0,
				  κ = 500.0,
				  σ = 50.0,
				  accept_r = 0,
				  accept_a = 0,
				  accept_κ = 0,
				  accept_K = 0,
				  accept_u0 = fill(0, N),
				  accept_σ = 0,
				  u = fill(0.0, m.T, m.N),
				  loglik = fill(0.0, N))

@btime sample_r!(pars, m)

log_K_star = m.μ_K + m.L_K * pars.α_K
p_star =  [pars.log_r, log_K_star, pars.a, pars.κ, m.λ]
@code_warntype process_all(p_star, pars.u0, m)

@code_warntype mcmc(m, pars, 1)

sum(chain["accept_K"])
sum(chain["accept_r"])

# Univariate plots

sum(chain["accept_sigma"][25001:end]) / 25000

plot(chain["sigma"])
histogram(chain["sigma"][25001:end], normalize = :true)
plot!(truncated(Normal(0, 100), 0.0, Inf))

sum(chain["accept_r"][25001:end]) / 25000

plot(chain["r"][1, :])
plot!(chain["r"][2, :])

density(exp.(chain["r"][1, 25001:end]), normalize = :true)
density!(chain["r"][2, 25001:end], normalize = :true)
plot!(Normal(log(0.2), 2.0))

sum(chain["accept_u0"], dims = 2) ./ 50000

plot(chain["u0"][end, 25001:end])

histogram(chain["u0"][1, 25001:end], normalize = :true)
plot!(Gamma(5, 200))

sum(chain["accept_K"][25001:end]) / 25000

plot(chain["K"][1, :])
plot!(chain["K"][2, :])
histogram(exp.(chain["K"][2, 25001:end]), normalize = :true)

sum(chain["accept_kappa"][25001:end]) / 25000

plot(chain["kappa"])
histogram(chain["kappa"][25001:end], normalize = :true)
plot!(Gamma(10, 50))

sum(chain["accept_a"][25001:end]) / 25000

plot(chain["a"])
histogram(chain["a"][25001:end], normalize = :true)
plot!(Gamma(1, 50))

# Dynamics
plot(1:26, chain["u"][:, 6, (end - 100):end], color = :gray, legend = false)
scatter!(1:26, m.z[:, 6])
hline!([exp(mean(chain["K"][1, :]))])

# Bivariate plots
scatter(chain["r"][25001:end], chain["K"][25001:end])
scatter(chain["r"][25001:end], chain["a"][25001:end])
scatter(chain["r"][25001:end], chain["kappa"][25001:end])

scatter(chain["K"][25001:end], chain["a"][25001:end])
scatter(chain["K"][225001:end], chain["kappa"][25001:end])

scatter(chain["a"][25001:end], chain["kappa"][25001:end])

scatter(chain["u0"][1, 25001:end], chain["K"][25001:end])
