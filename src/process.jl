# Code defining and solving the process model

# Rosenzweig-MacArthur model of single site prey dynamics
function prey!(du, u, p, t)
	@unpack log_r, log_K, a, κ, λ = p

	du[1] = exp(log_r[1]) * u[1] * (1.0 - u[1] / exp(log_K[1])) - a * u[1] * λ[1](t) / (u[1] + κ)
end

function process_site(p, u0, m)

	@unpack tspan, λ, N = m

	prob = ODEProblem(prey!, [u0], tspan, p)

	sol = solve(prob, Tsit5(), abstol = 1e-6, reltol = 1e-4, saveat=1.0)

	return sol[1, :]

end

# Rosenzweig-MacArthur model of multi-site prey dynamics
function prey_all!(du, u, p, t)
	@unpack r, a, κ, K, λ = p

	for i in 1:length(u)
		du[i] = r[i] * u[i] * (1.0 - u[i] / K) - a * λ[i](t) * u[i] / (u[i] + κ)
	end

end

# Solving ODE given parameters and λ
function process_all(p, u0, m)

	@unpack tspan = m

	prob = ODEProblem(prey_all!, u0, tspan, p)

	sol = solve(prob, Tsit5(), abstol = 1e-6, reltol = 1e-4, saveat=1.0)

	return sol[:, :]'

end
