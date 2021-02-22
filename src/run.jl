using DifferentialEquations
using Interpolations
using Parameters
using Distributions
using Plots
using StatsPlots
using RCall
using DataFrames

include("model.jl")
include("process.jl")
include("sample.jl")

# Reading in sea otter forcing time series
R"data <- readRDS('data/otter_array.rds')"
@rget data

# Reading in prey data
R"sag <- readRDS('data/SAG_array.rds')"
@rget sag

# Subsetting to sites with the prey observed
# obs = findall(sum(coalesce.(sag, 0.0), dims = 1)[1, :] .> 0)
data = data[:, 8]
sag = sag[:, 8]

# Number of sites
# N = size(data)[2]
N = 1

# Number of years
T = size(data)[1]

# Buidling interpolator to give the solver
λ = [interpolate(data[:, i], BSpline(Linear())) for i in 1:N]

# Simulating prey data
p = [0.2, 1000.0, 50.0, 500.0, λ]

prob = ODEProblem(prey!, rand(Normal(1000, 100), N), (1.0, 26.0), p)

sol = solve(prob, Tsit5(), saveat=1.0)

z = [rand(truncated(Normal(sol[i, t], 50), 0.0, Inf)) for t in 1:T, i in 1:N]

# Setting up model and parameter objects
m = model(z = sag,
          T = T,
		  tobs = 1:T,
          tspan = (1.0, 26.0),
		  N = N,
          λ = λ,
		  r_tune = 0.2,
		  r_prior = Gamma(1, 1),
		  a_tune = 20.0,
		  a_prior = Gamma(1, 50),
		  κ_tune = 100.0,
		  κ_prior = Gamma(50, 10),
		  K_tune = 100.0,
		  K_prior = Gamma(5, 100),
		  u0_tune = 500.0,
		  u0_prior = Gamma(5, 100)
		  )

pars = parameters(u0 = fill(600.0, N),
                  r = 0.2,
				  K = 600.0,
				  a = 100.0,
				  κ = 500.0,
				  accept_r = 0,
				  accept_a = 0,
				  accept_κ = 0,
				  accept_K = 0,
				  accept_u0 = 0,
				  u = sol,
				  loglik = 0.0)

chain = mcmc(m, pars, 50000)

# Univariate plots
sum(chain["accept_u0"][25001:end]) / 25000

plot(chain["u0"][1, :])
histogram(chain["u0"][1, 25001:end], normalize = :true)
plot!(Gamma(5, 100))

sum(chain["accept_K"][25001:end]) / 25000

plot(chain["K"])
histogram(chain["K"][25001:end], normalize = :true)
plot!(Gamma(5, 100))

sum(chain["accept_kappa"][25001:end]) / 25000

plot(chain["kappa"])
histogram(chain["kappa"][25001:end], normalize = :true)
plot!(Gamma(10, 50))

sum(chain["accept_a"][25001:end]) / 25000

plot(chain["a"])
histogram(chain["a"][25001:end], normalize = :true)
plot!(Gamma(1, 50))

sum(chain["accept_r"][25001:end]) / 25000

plot(chain["r"])
histogram(chain["r"][25001:end], normalize = :true)
plot!(Gamma(1, 1))

# Dynamics
plot(1:26, chain["u"][1, :, 49900:end], color = :gray, legend = false)
scatter!(1:26, m.z[:, 1])

# Bivariate plots
scatter(chain["r"][25001:end], chain["K"][25001:end])
scatter(chain["r"][25001:end], chain["a"][25001:end])
scatter(chain["r"][25001:end], chain["kappa"][25001:end])

scatter(chain["K"][25001:end], chain["a"][25001:end])
scatter(chain["K"][225001:end], chain["kappa"][25001:end])

scatter(chain["a"][25001:end], chain["kappa"][25001:end])

scatter(chain["u0"][1, 25001:end], chain["K"][25001:end])
