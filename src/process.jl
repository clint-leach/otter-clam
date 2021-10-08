# Rosenzweig-MacArthur model of multi-site prey dynamics
function prey_all!(du, u, p, t)
	@unpack r, a, κ, ν, λ = p

	for i in 1:length(u)
		du[i] = r[i] * u[i] - ν * u[i] ^ 2 - a[i] * λ[i](t) * u[i] ^ 2 / (u[i] ^ 2 + κ ^ 2)
	end
end

# Solving ODE given parameters and λ
function process_all(p, u0, m)

	@unpack tspan = m

	prob = ODEProblem(prey_all!, u0, tspan, p)

	sol = solve(prob, Tsit5(), abstol = 1e-8, reltol = 1e-6, saveat=1.0)

	return sol[:, :]'

end

# Solving ODE in prediction context (where state space is large)
function process_pred(p, u0, m)

	@unpack tspan = m

	u = fill(0.0, m.T, length(u0))

	site_p = DEparams([0.0], [0.0], 0.0, 0.0, [p.λ[1]])
	prob = ODEProblem(prey_all!, [0.0], tspan, site_p)

	for i in 1:length(u0)

		site_p = DEparams([p.r[i]], [p.a[i]], p.κ, p.ν, [p.λ[i]])

		prob = remake(prob; p = site_p, u0 = [u0[i]])

		sol = solve(prob, Tsit5(), abstol = 1e-8, reltol = 1e-6, saveat=1.0)

		u[:, i] = Array(sol)
	end

	return u

end
