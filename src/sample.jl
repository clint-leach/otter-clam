# Sampling scripts and helper functions

function likelihood(sol, m)

	@unpack z, tobs, N = m

	loglik = 0.0
	for i in 1:N
		for t in 1:length(tobs)
			if !ismissing(z[t, i])
				loglik += logpdf(truncated(Normal(sol[i, tobs[t]], 50.0), 0.0, Inf), z[t, i])
			end
		end
	end

	return loglik
end

# Sample prey population growth rate
function sample_r!(pars, m)

	@unpack r, accept_r, u0, K, a, κ, loglik, u = pars
	@unpack λ, r_tune, r_prior = m

	# Proposal
	forward_prop = truncated(Normal(r, r_tune), 0.0, Inf)
	r_star = rand(forward_prop)
	back_prop = truncated(Normal(r_star, r_tune), 0.0, Inf)

	# Proposal process model
	p_star =  [r_star, K, a, κ, λ]
	u_star = process(p_star, u0, m)

	# Proposal likelihood
	if u_star.t[end] < 26.0
		loglik_star = -Inf
	else
		loglik_star = likelihood(u_star, m)
	end

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(r_prior, r_star) + logpdf(back_prop, r)
	mh2 = loglik + logpdf(r_prior, r) + logpdf(forward_prop, r_star)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_r = 0
	else
		accept_r = 1
		r = r_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = r, accept_r, u, loglik
end

# Sample attack rate
function sample_a!(pars, m)

	@unpack r, accept_a, u0, K, a, κ, loglik, u = pars
	@unpack λ, a_tune, a_prior = m

	# Proposal
	forward_prop = truncated(Normal(a, a_tune), 0.0, Inf)
	a_star = rand(forward_prop)
	back_prop = truncated(Normal(a_star, a_tune), 0.0, Inf)

	# Proposal process model
	p_star =  [r, K, a_star, κ, λ]
	u_star = process(p_star, u0, m)

	# Proposal likelihood
	if u_star.t[end] < 26.0
		loglik_star = -Inf
	else
		loglik_star = likelihood(u_star, m)
	end

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(a_prior, a_star)  + logpdf(back_prop, a)
	mh2 = loglik + logpdf(a_prior, a) + logpdf(forward_prop, a_star)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_a = 0
	else
		accept_a = 1
		a = a_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = a, accept_a, u, loglik
end

# Sample functional response saturation constant
function sample_κ!(pars, m)

	@unpack r, accept_κ, u0, K, a, κ, loglik, u = pars
	@unpack λ, κ_tune, κ_prior = m

	# Proposal
	forward_prop = truncated(Normal(κ, κ_tune), 0.0, Inf)
	κ_star = rand(forward_prop)
	back_prop = truncated(Normal(κ_star, κ_tune), 0.0, Inf)

	# Proposal process model
	p_star =  [r, K, a, κ_star, λ]
	u_star = process(p_star, u0, m)

	# Proposal likelihood
	if u_star.t[end] < 26.0
		loglik_star = -Inf
	else
		loglik_star = likelihood(u_star, m)
	end

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(κ_prior, κ_star) + logpdf(back_prop, κ)
	mh2 = loglik + logpdf(κ_prior, κ) + logpdf(forward_prop, κ_star)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_κ = 0
	else
		accept_κ = 1
		κ = κ_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = κ, accept_κ, u, loglik
end

# Sample functional response saturation constant
function sample_K!(pars, m)

	@unpack r, accept_K, u0, K, a, κ, loglik, u, u0 = pars
	@unpack λ, K_tune, K_prior = m

	# Proposal
	forward_prop = truncated(Normal(K, K_tune), 0.0, Inf)
	K_star = rand(forward_prop)
	back_prop = truncated(Normal(K_star, K_tune), 0.0, Inf)

	# Proposal process model
	p_star =  [r, K_star, a, κ, λ]
	u_star = process(p_star, u0, m)

	# Proposal likelihood
	if u_star.t[end] < 26.0
		loglik_star = -Inf
	else
		loglik_star = likelihood(u_star, m)
	end

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(K_prior, K_star) + logpdf(back_prop, K)
	mh2 = loglik + logpdf(K_prior, K) + logpdf(forward_prop, K_star)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_K = 0
	else
		accept_K = 1
		K = K_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = K, accept_K, u, loglik
end

# Sample initial conditions
function sample_u0!(pars, m)

	@unpack r, accept_u0, u0, K, a, κ, loglik, u = pars
	@unpack λ, u0_tune, u0_prior, N = m

	forward_prop = truncated.(Normal.(u0, u0_tune), 0.0, Inf)
	u0_star = rand.(forward_prop)
	back_prop = truncated.(Normal.(u0_star, u0_tune), 0.0, Inf)

	# Proposal process model
	p =  [r, K, a, κ, λ]
	u_star = process(p, u0_star, m)

	# Proposal likelihood
	if u_star.t[end] < 26.0
		loglik_star = -Inf
	else
		loglik_star = likelihood(u_star, m)
	end

	# Computing the MH ratio
	mh1 = loglik_star + sum(logpdf.(u0_prior, u0_star)) + sum(logpdf.(back_prop, u0))
	mh2 = loglik + sum(logpdf.(u0_prior, u0)) +  sum(logpdf.(forward_prop, u0_star))

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_u0 = 0
	else
		accept_u0 = 1
		u0 = u0_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = u0, accept_u0, u, loglik
end

function mcmc(m, pars, nmcmc)

	chain = Dict("r" => fill(0.0, nmcmc),
	             "a" => fill(0.0, nmcmc),
				 "kappa" => fill(0.0, nmcmc),
				 "K" => fill(0.0, nmcmc),
				 "u0" => fill(0.0, N, nmcmc),
	             "accept_r" => fill(0, nmcmc),
				 "accept_a" => fill(0, nmcmc),
				 "accept_kappa" => fill(0, nmcmc),
				 "accept_K" => fill(0, nmcmc),
				 "accept_u0" => fill(0, nmcmc),
				 "u" => fill(0.0, m.N, m.T, nmcmc))

	# Initialize process and likelihood
	p =  [pars.r, pars.K, pars.a, pars.κ, m.λ]
	pars.u = process(p, pars.u0, m)
	pars.loglik = likelihood(pars.u, m)

	@progress for i in 1:nmcmc

		# Sampling

		sample_r!(pars, m)

		sample_a!(pars, m)

		sample_κ!(pars, m)

		sample_K!(pars, m)

		sample_u0!(pars, m)

		# Saving samples

		chain["r"][i] = pars.r
		chain["a"][i] = pars.a
		chain["kappa"][i] = pars.κ
		chain["K"][i] = pars.K
		chain["u0"][:, i] = pars.u0

		chain["accept_r"][i] = pars.accept_r
		chain["accept_a"][i] = pars.accept_a
		chain["accept_kappa"][i] = pars.accept_κ
		chain["accept_K"][i] = pars.accept_K
		chain["accept_u0"][i] = pars.accept_u0

		chain["u"][:, :, i] = [pars.u[j, t] for j in 1:m.N, t in 1:m.T]

	end

	return chain

end
