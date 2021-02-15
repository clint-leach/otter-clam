using DifferentialEquations
using Interpolations
using Parameters
using Distributions
using Plots
using StatsPlots
using RCall

include("model.jl")
include("process.jl")
include("sample.jl")

# Reading in sea otter forcing time series
R"data <- readRDS('data/otter_array.rds')"
@rget data

# Reading in prey data
R"sag <- readRDS('data/sag_array.rds')"
@rget sag

# Number of sites
N = size(data)[2]

# Number of years
T = size(data)[1]

# Buidling interpolator to give the solver
λ = [interpolate(data[:, i], BSpline(Linear())) for i in 1:N]

# Simulating prey data
p = [0.2, 1000.0, 50.0, 500.0, λ]

prob = ODEProblem(prey!, fill(1000, N), (1.0, 26.0), p)

sol = solve(prob, Tsit5(), saveat=1.0)

z = [rand(truncated(Normal(sol[i, t], 50), 0.0, Inf)) for t in 1:T, i in 1:N]

tobs = [5, 10, 15, 20, 25]
zobs = z[tobs, :]

# Setting up model and parameter objects
m = model(z = z,
          T = T,
		  tobs = 1:26,
          tspan = (1.0, 26.0),
		  N = N,
          λ = λ,
		  r_tune = 0.01,
		  r_prior = Gamma(1, 1),
		  a_tune = 0.1,
		  a_prior = Gamma(10, 10),
		  κ_tune = 2.0,
		  κ_prior = Gamma(5, 100),
		  K_tune = 0.1,
		  K_prior = Gamma(5, 100)
		  )

pars = parameters(u0 = 1000.0,
                  r = 0.2,
				  K = 1000.0,
				  a = 50.0,
				  κ = 500.0,
				  accept_r = 0,
				  accept_a = 0,
				  accept_κ = 0,
				  accept_K = 0,
				  u = sol,
				  loglik = 0.0)

chain = mcmc(m, pars, 10000)

sum(chain["accept_K"][5001:10000]) / 5000

plot(chain["K"])
histogram(chain["K"][5001:10000])

sum(chain["accept_kappa"][5001:10000]) / 5000

plot(chain["kappa"])
histogram(chain["kappa"][5001:10000])

sum(chain["accept_a"][5001:10000]) / 5000

plot(chain["a"])
histogram(chain["a"][5001:10000])

sum(chain["accept_r"][5001:10000]) / 5000

plot(chain["r"])
histogram(chain["r"][5001:10000])

scatter(1:26, z[:, 20])
# scatter!(tobs, zobs)
plot!(1:26, chain["u"][20, :, 9900:10000], color = :gray, legend = false)

scatter(chain["r"][5001:end], chain["K"][5001:end])
scatter(chain["r"][5001:end], chain["a"][5001:end])
scatter(chain["r"][5001:end], chain["kappa"][5001:end])

scatter(chain["K"][5001:end], chain["a"][5001:end])
scatter(chain["K"][5001:end], chain["kappa"][5001:end])

scatter(chain["a"][5001:end], chain["kappa"][5001:end])
