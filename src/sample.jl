# Sampling scripts and helper functions

function likelihood(u, m)

	@unpack z, T, N, nq = m

	loglik = fill(0.0, N)

	for i in 1:N
		for t in 1:T
			if !ismissing(z[t, i])
				loglik[i] += logpdf(Poisson(nq[i] * u[t, i]), z[t, i])
			end
		end
	end

	return loglik
end

# Sample z
function sample_z!(pars, m)

	@unpack T, N, nq = m
	@unpack z, u = pars

	for i in 1:N
		for t in 1:T
			z[t, i] = rand(Poisson(u[t, i]))
		end
	end

	@pack! pars = z
end

# Sample measurement variance
function sample_σ!(pars, m)

	@unpack loglik, u, σ, accept_σ, r = pars
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

	@unpack η_r, r, accept_r, u0, a, κ, loglik, u, K = pars
	@unpack λ, η_r_prior, r_tune, N, T = m

	# Proposal
	forward_prop = MvNormal(η_r, r_tune)
	η_r_star = rand(forward_prop)
	r_star = exp.(η_r_star)

	# Proposal process model
	p_star =  DEparams(r_star, a, κ, K, λ)
	u_star = process_all(p_star, u0, m)

	# Proposal likelihood
	if size(u_star, 1) < T
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, m)
	end

	# Computing the MH ratio
	mh1 = sum(loglik_star) + logpdf(η_r_prior, η_r_star)
	mh2 = sum(loglik) + logpdf(η_r_prior, η_r)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_r = 0
	else
		accept_r = 1
		η_r = η_r_star
		r = r_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = η_r, r, accept_r, u, u0, loglik
end

# Sample attack rate
function sample_a!(pars, m)

	@unpack r, accept_a, u0, a, κ, loglik, u, K = pars
	@unpack λ, a_tune, a_prior, N, T = m

	# Proposal
	forward_prop = truncated(Normal(a, a_tune), 0.0, Inf)
	a_star = rand(forward_prop)
	back_prop = truncated(Normal(a_star, a_tune), 0.0, Inf)

	# Proposal process model
	p_star =  DEparams(r, a_star, κ, K, λ)
	u_star = process_all(p_star, u0, m)

	# Proposal likelihood
	if size(u_star, 1) < T
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, m)
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

	@unpack r, accept_κ, u0, a, κ, loglik, u, K = pars
	@unpack λ, κ_tune, κ_prior, N, T = m

	# Proposal
	forward_prop = truncated(Normal(κ, κ_tune), 0.0, Inf)
	κ_star = rand(forward_prop)
	back_prop = truncated(Normal(κ_star, κ_tune), 0.0, Inf)

	# Proposal process model
	p_star =  DEparams(r, a, κ_star, K, λ)
	u_star = process_all(p_star, u0, m)

	# Proposal likelihood
	if size(u_star, 1) < T
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, m)
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

function sample_K!(pars, m)

	@unpack r, K, accept_K, u0, a, κ, loglik, u, η_0 = pars
	@unpack λ, K_tune, K_prior, N, T = m

	# Proposal
	forward_prop = truncated(Normal(K, K_tune), 0.0, Inf)
	K_star = rand(forward_prop)
	back_prop = truncated(Normal(K_star, K_tune), 0.0, Inf)

	# Proposal process model
	p_star =  DEparams(r, a, κ, K_star, λ)
	u0_star = K_star * logistic.(η_0)
	u_star = process_all(p_star, u0_star, m)

	# Proposal likelihood
	if size(u_star, 1) < T
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, m)
	end

	# Computing the MH ratio
	mh1 = sum(loglik_star) + logpdf(K_prior, K_star) + logpdf(back_prop, K)
	mh2 = sum(loglik) + logpdf(K_prior, K) + logpdf(forward_prop, K_star)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_K = 0
	else
		accept_K = 1
		K = K_star
		u = u_star
		u0 = u0_star
		loglik = loglik_star
	end

	@pack! pars = K, accept_K, u, u0, loglik
end

# Sample functional response saturation constant
function sample_u0!(pars, m)

	@unpack r, accept_u0, u0, η_0, a, κ, loglik, u, K = pars
	@unpack λ, η_0_prior, u0_tune, N, T = m

	# Proposal
	forward_prop = MvNormal(η_0, u0_tune)
	η_0_star = rand(forward_prop)
	u0_star = K * logistic.(η_0_star)

	# Proposal process model
	p =  DEparams(r, a, κ, K, λ)
	u_star = process_all(p, u0_star, m)

	# Proposal likelihood
	if size(u_star, 1) < T
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, m)
	end

	# Computing the MH ratio
	mh1 = sum(loglik_star) + logpdf(η_0_prior, η_0_star)
	mh2 = sum(loglik) + logpdf(η_0_prior, η_0)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_u0 = 0
	else
		accept_u0 = 1
		η_0 = η_0_star
		u0 = u0_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = η_0, accept_u0, u0, u, loglik
end


# Conjugate Gibbs updates of regression coefficients on r
function sample_β_r!(pars, m)

	@unpack η_r = pars
	@unpack X, Ω_r, Ω_β_r, μ_β_r = m

	# Sample β_r
	A = Symmetric(X' * Ω_r * X + Ω_β_r)
	A_inv = inv(A)

	b = X' * Ω_r * η_r + Ω_β_r * μ_β_r

	β_r = rand(MvNormal(A_inv * b, A_inv))

	@pack! pars = β_r

end

# Conjugate Gibbs updates of regression coefficients on K
function sample_β_0!(pars, m)

	@unpack η_0 = pars
	@unpack X, Ω_0, Ω_β_0, μ_β_0 = m

	# Sample β_r
	A = Symmetric(X' * Ω_0 * X + Ω_β_0)
	A_inv = inv(A)

	b = X' * Ω_0 * η_0 + Ω_β_0 * μ_β_0

	β_0 = rand(MvNormal(A_inv * b, A_inv))

	@pack! pars = β_0

end

function mcmc(m, pars, keep_every, nburn, nmcmc)

	nkeep = floor(Int64, nmcmc / keep_every)

	chain = Dict("r" => fill(0.0, m.N, nkeep),
				 "eta_r"=> fill(0.0, m.N, nkeep),
	             "a" => fill(0.0, nkeep),
				 "K" => fill(0.0, nkeep),
				 "kappa" => fill(0.0, nkeep),
				 "eta_0" => fill(0.0, m.N, nkeep),
				 "u0" => fill(0.0, m.N, nkeep),
				 "beta_r" => fill(0.0, m.p, nkeep),
				 "beta_0" => fill(0.0, m.p, nkeep),
	             "accept_r" => fill(0, nkeep),
				 "accept_a" => fill(0, nkeep),
				 "accept_kappa" => fill(0, nkeep),
				 "accept_K" => fill(0, nkeep),
				 "accept_u0" => fill(0, nkeep),
				 "u" => fill(0.0, m.T, m.N, nkeep),
				 "zpred" => fill(0.0, m.T, m.N, nkeep))

	# Initialize process and likelihood
	p = DEparams(pars.r, pars.a, pars.κ, pars.K, m.λ)
	pars.u = process_all(p, pars.u0, m)
	pars.loglik = likelihood(pars.u, m)

	# Burn-in
	for i in 1:nburn

		if i % 500 == 0
			@show i
		end

		# Sampling

		sample_K!(pars, m)

		sample_a!(pars, m)

		sample_κ!(pars, m)

		sample_u0!(pars, m)

		sample_r!(pars, m)

	end

	for i in 1:nmcmc

		if i % 500 == 0
			@show i
		end

		# Sampling
		sample_K!(pars, m)

		sample_a!(pars, m)

		sample_κ!(pars, m)

		sample_u0!(pars, m)

		sample_β_0!(pars, m)

		sample_r!(pars, m)

		sample_β_r!(pars, m)

		sample_z!(pars, m)

		# Saving samples
		if i % keep_every == 0
			idx = floor(Int64, i / keep_every)

			chain["r"][:, idx] = pars.r
			chain["eta_r"][:, idx] = pars.η_r
			chain["a"][idx] = pars.a
			chain["K"][idx] = pars.K
			chain["kappa"][idx] = pars.κ
			chain["eta_0"][:, idx] = pars.η_0
			chain["u0"][:, idx] = pars.u0
			chain["beta_r"][:, idx] = pars.β_r
			chain["beta_0"][:, idx] = pars.β_0

			chain["accept_r"][idx] = pars.accept_r
			chain["accept_a"][idx] = pars.accept_a
			chain["accept_kappa"][idx] = pars.accept_κ
			chain["accept_K"][idx] = pars.accept_K
			chain["accept_u0"][idx] = pars.accept_u0

			chain["u"][:, :, idx] = pars.u
			chain["zpred"][:, :, idx] = pars.z
		end
	end

	return chain

end
