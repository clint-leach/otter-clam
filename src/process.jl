# Code defining and solving the process model

# Rosenzweig-MacArthur model of prey dynamics
function prey!(du, u, p, t)
	log_r, log_K, a, κ, λ = p

	for i in 1:length(u)
		du[i] = exp(log_r[i]) * u[i] * (1.0 - u[i] / exp(log_K[i])) - a * u[i] * λ[i](t) / (u[i] + κ)
	end
end

# Solving ODE given parameters and λ
function process(p, u0, m)

	@unpack tspan, λ, N = m

	prob = ODEProblem(prey!, u0, tspan, p)

	sol = solve(prob, Tsit5(), saveat=1.0)

	return sol

end
