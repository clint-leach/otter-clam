# Sampling scripts and helper functions

function likelihood(u, m)

	loglik = 0.0
	for t in 1:m.T
		loglik += logpdf(truncated(Normal(u[t], 50.0), 0.0, Inf), z[t])
	end

	return loglik
end

# Sample prey population growth rate
function sample_r!(pars, m)

	@unpack r, accept_r, u0, K, a, κ, loglik, u = pars
	@unpack λ, r_tune, r_prior = m

	# Proposal
	r_star = rand(Gamma(r / r_tune, r_tune))

	# Proposal process model
	p_star =  [r_star, K, a, κ, λ]
	u_star = process(p_star, u0, m)

	# Proposal likelihood
	loglik_star = likelihood(u_star, m)

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(r_prior, r_star) + logpdf(Gamma(r_star / r_tune, r_tune), r)
	mh2 = loglik + logpdf(r_prior, r) +  logpdf(Gamma(r / r_tune, r_tune), r_star)

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
	a_star = rand(Gamma(a / a_tune, a_tune))

	# Proposal process model
	p_star =  [r, K, a_star, κ, λ]
	u_star = process(p_star, u0, m)

	# Proposal likelihood
	loglik_star = likelihood(u_star, m)

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(a_prior, a_star) + logpdf(Gamma(a_star / a_tune, a_tune), a)
	mh2 = loglik + logpdf(a_prior, a) +  logpdf(Gamma(a / a_tune, a_tune), a_star)

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
	κ_star = rand(Gamma(κ / κ_tune, κ_tune))

	# Proposal process model
	p_star =  [r, K, a, κ_star, λ]
	u_star = process(p_star, u0, m)

	# Proposal likelihood
	loglik_star = likelihood(u_star, m)

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(κ_prior, κ_star) + logpdf(Gamma(κ_star / κ_tune, κ_tune), κ)
	mh2 = loglik + logpdf(κ_prior, κ) +  logpdf(Gamma(κ / κ_tune, κ_tune), κ_star)

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

	@unpack r, accept_K, u0, K, a, κ, loglik, u = pars
	@unpack λ, K_tune, K_prior = m

	# Proposal
	K_star = rand(Gamma(K / K_tune, K_tune))

	# Proposal process model
	p_star =  [r, K_star, a, κ, λ]
	u_star = process(p_star, u0, m)

	# Proposal likelihood
	loglik_star = likelihood(u_star, m)

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(K_prior, K_star) + logpdf(Gamma(K_star / K_tune, K_tune), K)
	mh2 = loglik + logpdf(K_prior, K) +  logpdf(Gamma(K / K_tune, K_tune), K_star)

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

function mcmc(m, pars, nmcmc)

	chain = Dict("r" => fill(0.0, nmcmc),
	             "a" => fill(0.0, nmcmc),
				 "kappa" => fill(0.0, nmcmc),
				 "K" => fill(0.0, nmcmc),
	             "accept_r" => fill(0, nmcmc),
				 "accept_a" => fill(0, nmcmc),
				 "accept_kappa" => fill(0, nmcmc),
				 "accept_K" => fill(0, nmcmc),
				 "u" => fill(0.0, 26, nmcmc))

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

		# Saving samples

		chain["r"][i] = pars.r
		chain["a"][i] = pars.a
		chain["kappa"][i] = pars.κ
		chain["K"][i] = pars.K

		chain["accept_r"][i] = pars.accept_r
		chain["accept_a"][i] = pars.accept_a
		chain["accept_kappa"][i] = pars.accept_κ
		chain["accept_K"][i] = pars.accept_K

		chain["u"][:, i] = pars.u

	end

	return chain

end
