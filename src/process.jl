# Code defining and solving the process model

# Rosenzweig-MacArthur model of prey dynamics
function prey!(du, u, p, t)
	r, K, a, κ, λ = p

	du[1] = r * u[1] * (1.0 - u[1] / K) - a * u[1] * λ(t) / (u[1] + κ)
end

# Solving ODE given parameters and λ
function process(p, u0, m)

	@unpack tspan, λ = m

	prob = ODEProblem(prey!, [u0], tspan, p)

	sol = solve(prob, Tsit5(), saveat=1.0)

	return sol[1, :]

end
