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

# Setting up model and parameter objects
m = model(z = z,
          T = length(z),
          tspan = (1.0, 26.0),
          λ = λ,
		  r_tune = 0.001,
		  r_prior = Gamma(1, 1)
		  )

pars = parameters(u0 = 1000.0,
                  r = 0.5,
				  K = 1000.0,
				  a = 50.0,
				  κ = 500.0,
				  accept_r = 0,
				  u = fill(0.0, 26),
				  loglik = 0.0)

sample_r!(pars, m)

chain = mcmc(m, pars, 1000)

sum(chain["accept_r"])

plot(chain["r"])
histogram(chain["r"][501:1000])

scatter(1:26, z)
plot!(1:26, chain["u"][:, end])
