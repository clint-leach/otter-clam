# Rosenzweig-MacArthur model of multi-site prey dynamics
function prey_all!(du, u, p, t)
	@unpack r, a, κ, K, λ = p

	for i in 1:length(u)
		du[i] = r[i] * u[i] * (1.0 - u[i] / K) - a * λ[i](t) * u[i] ^ 2 / (u[i] ^ 2 + κ ^ 2)
	end
end

# Solving ODE given parameters and λ
function process_all(p, u0, m)

	@unpack tspan = m

	prob = ODEProblem(prey_all!, u0, tspan, p)

	sol = solve(prob, AutoTsit5(Rosenbrock23()), abstol = 1e-4, reltol = 1e-4, saveat=1.0)

	return sol[:, :]'

end

# Solving ODE in prediction context (where state space is large)
function process_pred(p, u0, m)

	@unpack tspan = m

	u = fill(0.0, m.T, length(u0))

	for i in 1:length(u0)

		site_p = DEparams([p.r[i]], p.a, p.κ, p.K, [p.λ[i]])

		prob = ODEProblem(prey_all!, [u0[i]], tspan, site_p)

		sol = solve(prob, AutoVern7(Rosenbrock23()), abstol = 1e-10, reltol = 1e-4, saveat=1.0)

		u[:, i] = Array(sol)
	end

	return u

end
