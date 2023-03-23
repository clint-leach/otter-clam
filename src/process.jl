# Rosenzweig-MacArthur model of multi-site prey dynamics
function prey_all!(du, u, p, t)
	@unpack r, a, ν, λ = p

	for i in 1:length(u)
		du[i] = ν * (r[i] - u[i]) - a[i] * λ[i](t) * u[i]
	end
end

# Solving ODE given parameters and λ
function process_all(p, u0, m)

	@unpack tspan = m

	prob = ODEProblem(prey_all!, u0, tspan, p)

	sol = solve(prob, Tsit5(), abstol = 1e-8, reltol = 1e-6, saveat=1.0)

	return sol[:, :]'

end
