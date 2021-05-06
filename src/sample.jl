# Sampling scripts and helper functions

function likelihood(u, σ, m)

	@unpack z, T, N = m

	loglik = fill(0.0, N)

	for i in 1:N
		for t in 1:T
			if !ismissing(z[t, i])
				loglik[i] += logpdf(truncated(Normal(exp(u[t, i]), σ), 0.0, Inf), z[t, i])
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

	@unpack log_r, accept_r, u0, log_K, a, κ, loglik, u, σ = pars
	@unpack λ, log_r_prior, r_tune = m

	# Proposal
	forward_prop = MvNormal(log_r, r_tune)
	log_r_star = rand(forward_prop)

	# Proposal process model
	p_star =  DEparams(log_r_star, log_K, a, κ, λ)
	u_star = process_all(p_star, u0, m)

	# Proposal likelihood
	if size(u_star, 1) < 26
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, σ, m)
	end

	# Computing the MH ratio
	mh1 = sum(loglik_star) + logpdf(log_r_prior, log_r_star)
	mh2 = sum(loglik) + logpdf(log_r_prior, log_r)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_r = 0
	else
		accept_r = 1
		log_r = log_r_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = log_r, accept_r, u, loglik
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

	@unpack log_r, accept_K, u0, log_K, a, κ, loglik, u, σ = pars
	@unpack λ, log_K_prior, K_tune = m

	# Proposal
	forward_prop = MvNormal(log_K, K_tune)
	log_K_star = rand(forward_prop)

	# Proposal process model
	u0_star = log_K_star
	p_star =  DEparams(log_r, log_K_star, a, κ, λ)
	u_star = process_all(p_star, u0_star, m)

	# Proposal likelihood
	if size(u_star, 1) < 26
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, σ, m)
	end

	# Computing the MH ratio
	mh1 = sum(loglik_star) + logpdf(log_K_prior, log_K_star)
	mh2 = sum(loglik) + logpdf(log_K_prior, log_K)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_K = 0
	else
		accept_K = 1
		log_K = log_K_star
		u0 = u0_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = log_K, accept_K, u0, u, loglik
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

# Sample initial condition regression coefficients
function sample_β!(pars, m)

	@unpack log_r, accept_β, β, log_K, a, κ, loglik, u, σ, u0 = pars
	@unpack λ, β_tune, β_prior, N, X = m

	# Proposal
	forward_prop = MvNormal(β, β_tune)
	β_star = rand(forward_prop)

	# Proposal process model
	u0_star = exp.(X * β_star)
	p = DEparams(log_r, log_K, a, κ, λ)
	u_star = process_all(p, u0_star, m)

	# Proposal likelihood
	if size(u_star, 1) < 26
		loglik_star = fill(-Inf, N)
	else
		loglik_star = likelihood(u_star, σ, m)
	end

	# Computing the MH ratio
	mh1 = sum(loglik_star) + logpdf(β_prior, β_star)
	mh2 = sum(loglik) + logpdf(β_prior, β)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_β = 0
	else
		accept_β = 1
		β = β_star
		u = u_star
		loglik = loglik_star
		u0 = u0_star
	end


	@pack! pars = u0, accept_β, β, u, loglik
end

# Conjugate Gibbs updates of regression coefficients on r
function sample_β_r!(pars, m)

	@unpack log_r = pars
	@unpack X, Ω_r, Ω_β_r, μ_β_r = m

	# Sample β_r
	A = Symmetric(X' * Ω_r * X + Ω_β_r)
	A_inv = inv(A)

	b = X' * Ω_r * log_r + Ω_β_r * μ_β_r

	β_r = rand(MvNormal(A_inv * b, A_inv))

	@pack! pars = β_r

end

# Conjugate Gibbs updates of regression coefficients on K
function sample_β_K!(pars, m)

	@unpack log_K = pars
	@unpack X, Ω_K, Ω_β_K, μ_β_K = m

	# Sample β_r
	A = Symmetric(X' * Ω_K * X + Ω_β_K)
	A_inv = inv(A)

	b = X' * Ω_K * log_K + Ω_β_K * μ_β_K

	β_K = rand(MvNormal(A_inv * b, A_inv))

	@pack! pars = β_K

end

function mcmc(m, pars, nmcmc)

	chain = Dict("r" => fill(0.0, m.N, nmcmc),
	             "a" => fill(0.0, nmcmc),
				 "kappa" => fill(0.0, nmcmc),
				 "K" => fill(0.0, m.N, nmcmc),
				 "u0" => fill(0.0, m.N, nmcmc),
				 "sigma" => fill(0.0, nmcmc),
				 "beta_r" => fill(0.0, m.p, nmcmc),
				 "beta_K" => fill(0.0, m.p, nmcmc),
	             "accept_r" => fill(0, nmcmc),
				 "accept_a" => fill(0, nmcmc),
				 "accept_kappa" => fill(0, nmcmc),
				 "accept_K" => fill(0, nmcmc),
				 "accept_beta" => fill(0, nmcmc),
				 "accept_u0" => fill(0, m.N, nmcmc),
				 "accept_sigma" => fill(0, nmcmc),
				 "u" => fill(0.0, m.T, m.N, nmcmc))

	# Initialize process and likelihood
	pars.u0 = pars.log_K
	p =  DEparams(pars.log_r, pars.log_K, pars.a, pars.κ, m.λ)
	pars.u = process_all(p, pars.u0, m)
	pars.loglik = likelihood(pars.u, pars.σ, m)

	for i in 1:nmcmc

		if i % 100 == 0
			@show i
		end

		# Sampling

		sample_a!(pars, m)

		# sample_κ!(pars, m)

		sample_K!(pars, m)

		sample_β_K!(pars, m)

		sample_r!(pars, m)

		sample_β_r!(pars, m)

		sample_σ!(pars, m)

		# Saving samples

		chain["r"][:, i] = pars.log_r
		chain["a"][i] = pars.a
		chain["kappa"][i] = pars.κ
		chain["K"][:, i] = pars.log_K
		chain["u0"][:, i] = pars.u0
		chain["sigma"][i] = pars.σ
		chain["beta_r"][:, i] = pars.β_r
		chain["beta_K"][:, i] = pars.β_K

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
