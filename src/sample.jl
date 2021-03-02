# Sampling scripts and helper functions

function likelihood(u, σ, m)

	@unpack z, tobs, N = m

	loglik = fill(0.0, N)
	for i in 1:N
		for t in 1:length(tobs)
			if !ismissing(z[t, i])
				loglik[i] += logpdf(truncated(Normal(u[tobs[t], i], σ), 0.0, Inf), z[tobs[t], i])
			end
		end
	end

	return loglik
end

function likelihood(u, σ, m, i)

	@unpack z, tobs, N = m

	loglik = 0.0

	for t in 1:length(tobs)
		if !ismissing(z[t, i])
			loglik += logpdf(truncated(Normal(u[tobs[t]], σ), 0.0, Inf), z[tobs[t], i])
		end
	end

	return loglik
end

# Sample measurement variance
function sample_σ!(pars, m)

	@unpack loglik, u, σ, accept_σ = pars
	@unpack σ_prior, σ_tune = m

	# Proposal
	forward_prop = truncated(Normal(σ, σ_tune), 0.0, Inf)
	σ_star = rand(forward_prop)
	back_prop = truncated(Normal(σ_star, σ_tune), 0.0, Inf)

	# Proposal likelihood
	loglik_star = likelihood(u, σ_star, m)

	# Computing the MH ratio
	mh1 = sum(loglik_star) + logpdf(σ_prior, σ_star)  + logpdf(back_prop, σ)
	mh2 = sum(loglik) + logpdf(σ_prior, σ) + logpdf(forward_prop, σ_star)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_σ = 0
	else
		accept_σ = 1
		σ = σ_star
		loglik = loglik_star
	end

	@pack! pars = σ, accept_σ, loglik

end

# Sample prey population growth rate
function sample_r!(pars, m)

	@unpack log_r, α_r, accept_r, u0, log_K, a, κ, loglik, u, σ = pars
	@unpack λ, α_r_tune, α_r_prior, L_r, μ_r = m

	# Proposal
	forward_prop = MvNormal(α_r, α_r_tune)
	α_r_star = rand(forward_prop)

	# Proposal process model
	log_r_star = μ_r + L_r * α_r_star
	p_star =  DEparams(log_r_star, log_K, a, κ, λ)
	u_star = process_all(p_star, u0, m)

	# Proposal likelihood
	if size(u_star, 1) < 26
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, σ, m)
	end

	# Computing the MH ratio
	mh1 = sum(loglik_star) + logpdf(α_r_prior, α_r_star)
	mh2 = sum(loglik) + logpdf(α_r_prior, α_r)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_r = 0
	else
		accept_r = 1
		log_r = log_r_star
		α_r = α_r_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = log_r, accept_r, α_r, u, loglik
end

# Sample attack rate
function sample_a!(pars, m)

	@unpack log_r, accept_a, u0, log_K, a, κ, loglik, u, σ = pars
	@unpack λ, a_tune, a_prior = m

	# Proposal
	forward_prop = truncated(Normal(a, a_tune), 0.0, Inf)
	a_star = rand(forward_prop)
	back_prop = truncated(Normal(a_star, a_tune), 0.0, Inf)

	# Proposal process model
	p_star =  DEparams(log_r, log_K, a_star, κ, λ)
	u_star = process_all(p_star, u0, m)

	# Proposal likelihood
	if size(u_star, 1) < 26
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, σ, m)
	end

	# Computing the MH ratio
	mh1 = sum(loglik_star) + logpdf(a_prior, a_star)  + logpdf(back_prop, a)
	mh2 = sum(loglik) + logpdf(a_prior, a) + logpdf(forward_prop, a_star)

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

	@unpack log_r, accept_κ, u0, log_K, a, κ, loglik, u, σ = pars
	@unpack λ, κ_tune, κ_prior = m

	# Proposal
	forward_prop = truncated(Normal(κ, κ_tune), 0.0, Inf)
	κ_star = rand(forward_prop)
	back_prop = truncated(Normal(κ_star, κ_tune), 0.0, Inf)

	# Proposal process model
	p_star =  DEParams(log_r, log_K, a, κ_star, λ)
	u_star = process_all(p_star, u0, m)

	# Proposal likelihood
	if size(ustar, 1) < 26
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, σ, m)
	end

	# Computing the MH ratio
	mh1 = sum(loglik_star) + logpdf(κ_prior, κ_star) + logpdf(back_prop, κ)
	mh2 = sum(loglik) + logpdf(κ_prior, κ) + logpdf(forward_prop, κ_star)

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

	@unpack log_r, accept_K, u0, log_K, α_K, a, κ, loglik, u, σ = pars
	@unpack λ, α_K_tune, α_K_prior, μ_K, L_K = m

	# Proposal
	forward_prop = MvNormal(α_K, α_K_tune)
	α_K_star = rand(forward_prop)

	# Proposal process model
	log_K_star = μ_K + L_K * α_K_star
	p_star =  DEparams(log_r, log_K_star, a, κ, λ)
	u_star = process_all(p_star, u0, m)

	# Proposal likelihood
	if size(u_star, 1) < 26
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, σ, m)
	end

	# Computing the MH ratio
	mh1 = sum(loglik_star) + logpdf(α_K_prior, α_K_star)
	mh2 = sum(loglik) + logpdf(α_K_prior, α_K)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_K = 0
	else
		accept_K = 1
		log_K = log_K_star
		α_K = α_K_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = log_K, accept_K, α_K, u, loglik
end

# Sample initial conditions
function sample_u0!(pars, m)

	@unpack log_r, accept_u0, u0, log_K, a, κ, loglik, u, σ = pars
	@unpack λ, u0_tune, u0_prior, N = m

	for i in 1:N

		forward_prop = truncated(Normal(u0[i], u0_tune), 0.0, Inf)
		u0_star = rand(forward_prop)
		back_prop = truncated(Normal(u0_star, u0_tune), 0.0, Inf)

		# Proposal process model
		p =  DEparams(log_r[i, :], log_K[i, :], a, κ, λ[i, :])
		u_star = process_site(p, u0_star, m)

		# Proposal likelihood
		if length(u_star) < 26.0
			loglik_star = -Inf
		else
			loglik_star = likelihood(u_star, σ, m, i)
		end

		# Computing the MH ratio
		mh1 = loglik_star + sum(logpdf(u0_prior, u0_star)) + sum(logpdf(back_prop, u0[i]))
		mh2 = loglik[i] + sum(logpdf(u0_prior, u0[i])) +  sum(logpdf(forward_prop, u0_star))

		# Accept/reject
		prob = exp(mh1 - mh2)
		if rand() > prob
			accept_u0[i] = 0
		else
			accept_u0[i] = 1
			u0[i] = u0_star
			u[:, i] = u_star
			loglik[i] = loglik_star
		end
	end

	@pack! pars = u0, accept_u0, u, loglik
end

function mcmc(m, pars, nmcmc)

	chain = Dict("r" => fill(0.0, m.N, nmcmc),
	             "a" => fill(0.0, nmcmc),
				 "kappa" => fill(0.0, nmcmc),
				 "K" => fill(0.0, m.N, nmcmc),
				 "u0" => fill(0.0, m.N, nmcmc),
				 "sigma" => fill(0.0, nmcmc),
	             "accept_r" => fill(0, nmcmc),
				 "accept_a" => fill(0, nmcmc),
				 "accept_kappa" => fill(0, nmcmc),
				 "accept_K" => fill(0, nmcmc),
				 "accept_u0" => fill(0, m.N, nmcmc),
				 "accept_sigma" => fill(0, nmcmc),
				 "u" => fill(0.0, m.T, m.N, nmcmc))

	# Initialize process and likelihood
	p =  DEparams(pars.log_r, pars.log_K, pars.a, pars.κ, m.λ)
	pars.u = process_all(p, pars.u0, m)
	pars.loglik = likelihood(pars.u, pars.σ, m)

	@progress for i in 1:nmcmc

		# Sampling

		sample_r!(pars, m)

		sample_a!(pars, m)

		# sample_κ!(pars, m)

		sample_K!(pars, m)

		sample_u0!(pars, m)

		sample_σ!(pars, m)

		# Saving samples

		chain["r"][:, i] = pars.log_r
		chain["a"][i] = pars.a
		chain["kappa"][i] = pars.κ
		chain["K"][:, i] = pars.log_K
		chain["u0"][:, i] = pars.u0
		chain["sigma"][i] = pars.σ

		chain["accept_r"][i] = pars.accept_r
		chain["accept_a"][i] = pars.accept_a
		chain["accept_kappa"][i] = pars.accept_κ
		chain["accept_K"][i] = pars.accept_K
		chain["accept_u0"][:, i] = pars.accept_u0
		chain["accept_sigma"][i] = pars.accept_σ

		chain["u"][:, :, i] = pars.u

	end

	return chain

end
