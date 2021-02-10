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

# Reading in real data and forcing the prey model
R"data <- readRDS('data/otters.rds')"
@rget data

# Filtering down to a single site
site = filter(row -> row.site == "221", data)

# Buidling interpolator to give the solver
lambda = Array{Float64}(site[!, :lambda])
λ = interpolate(lambda, BSpline(Linear()))

# Simulating prey data
p = [0.2, 1000.0, 50.0, 500.0, λ]

prob = ODEProblem(prey!, [1000.0], (1.0, 26.0), p)

sol = solve(prob, Tsit5(), saveat=1.0)

z = rand.(truncated.(Normal.(sol[1, :], 50), 0.0, Inf))

tobs = [5, 10, 15, 20, 25]
zobs = z[tobs]

# Setting up model and parameter objects
m = model(z = zobs,
          T = length(λ),
		  tobs = tobs,
          tspan = (1.0, 26.0),
          λ = λ,
		  r_tune = 0.1,
		  r_prior = Gamma(1, 1),
		  a_tune = 0.2,
		  a_prior = Gamma(10, 10),
		  κ_tune = 10.0,
		  κ_prior = Gamma(5, 100),
		  K_tune = 5.0,
		  K_prior = Gamma(5, 100)
		  )

pars = parameters(u0 = 1000.0,
                  r = 0.2,
				  K = 500.0,
				  a = 50.0,
				  κ = 500.0,
				  accept_r = 0,
				  accept_a = 0,
				  accept_κ = 0,
				  accept_K = 0,
				  u = fill(0.0, 26),
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

scatter(1:26, z)
scatter!(tobs, zobs)
plot!(1:26, chain["u"][:, 9900:10000], color = :gray, legend = false)

scatter(chain["r"][5001:end], chain["K"][5001:end])
scatter(chain["r"][5001:end], chain["a"][5001:end])
scatter(chain["r"][5001:end], chain["kappa"][5001:end])

scatter(chain["K"][5001:end], chain["a"][5001:end])
scatter(chain["K"][5001:end], chain["kappa"][5001:end])

scatter(chain["a"][5001:end], chain["kappa"][5001:end])
