# Sampling scripts and helper functions

function likelihood(u, σ, m)

	@unpack z, T, N, nq = m

	loglik = 0.0

	for i in 1:N
		for t in 1:T
			if !ismissing(z[1, t, i])
				n = σ / (1 - σ) * u[t, i]
				loglik += sum(logpdf.(NegativeBinomial(n, σ), z[1:nq[i], t, i]))
			end
		end
	end

	return loglik
end

# Sample z
function sample_z!(pars, m)

	@unpack T, N = m
	@unpack z, u, σ = pars

	for i in 1:N
		for t in 1:T
			n = σ / (1 - σ) * u[t, i]
			z[t, i] = rand(NegativeBinomial(n, σ))
		end
	end

	@pack! pars = z
end

# Sample measurement variance
function sample_σ!(pars, m)

	@unpack loglik, u, σ, accept_σ, r = pars
	@unpack σ_prior, σ_tune = m

	# Proposal
	forward_prop = truncated(Normal(σ, σ_tune), 0.0, 1.0)
	σ_star = rand(forward_prop)
	back_prop = truncated(Normal(σ_star, σ_tune), 0.0, 1.0)

	# Proposal likelihood
	loglik_star = likelihood(u, σ_star, m)

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(σ_prior, σ_star)  + logpdf(back_prop, σ)
	mh2 = loglik + logpdf(σ_prior, σ) + logpdf(forward_prop, σ_star)

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

	@unpack r, η_r, accept_r, a, κ, loglik, u, ν, σ, u0, η_0 = pars
	@unpack λ, η_r_prior, r_tune, N, T = m

	# Proposal
	forward_prop = MvNormal(η_r, r_tune)
	η_r_star = rand(forward_prop)
	r_star = exp.(η_r_star)
	u0_star = logistic.(η_0) .* r_star ./ ν

	# Proposal process model
	p_star =  DEparams(r_star, a, κ, ν, λ)
	u_star = process_all(p_star, u0_star, m)

	# Proposal likelihood
	if size(u_star, 1) < T
		loglik_star = -Inf
	else
		loglik_star = likelihood(u_star, σ, m)
	end

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(η_r_prior, η_r_star)
	mh2 = loglik + logpdf(η_r_prior, η_r)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_r = 0
	else
		accept_r = 1
		η_r = η_r_star
		r = r_star
		u = u_star
		u0 = u0_star
		loglik = loglik_star
	end

	@pack! pars = η_r, r, accept_r, u, u0, loglik
end

# Sample attack rate
function sample_a!(pars, m)

	@unpack r, accept_a, u0, a, η_a, κ, loglik, u, ν, σ = pars
	@unpack λ, a_tune, η_a_prior, N, T = m

	# Proposal
	forward_prop = MvNormal(η_a, a_tune)
	η_a_star = rand(forward_prop)
	a_star = exp.(η_a_star)

	# Proposal process model
	p_star =  DEparams(r, a_star, κ, ν, λ)
	u_star = process_all(p_star, u0, m)

	# Proposal likelihood
	if size(u_star, 1) < T
		loglik_star = -Inf
	else
		loglik_star = likelihood(u_star, σ, m)
	end

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(η_a_prior, η_a_star)
	mh2 = loglik + logpdf(η_a_prior, η_a)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_a = 0
	else
		accept_a = 1
		a = a_star
		η_a = η_a_star
		u = u_star
		loglik = loglik_star
	end

	@pack! pars = a, η_a, accept_a, u, loglik
end

# Sample functional response saturation constant
function sample_κ!(pars, m)

	@unpack r, accept_κ, u0, a, κ, loglik, u, ν, σ = pars
	@unpack λ, κ_tune, κ_prior, N, T = m

	# Proposal
	forward_prop = Gamma(κ / κ_tune, κ_tune)
	κ_star = rand(forward_prop)
	back_prop = Gamma(κ_star / κ_tune, κ_tune)

	# Proposal process model
	p_star =  DEparams(r, a, κ_star, ν, λ)
	u_star = process_all(p_star, u0, m)

	# Proposal likelihood
	if size(u_star, 1) < T
		loglik_star = -Inf
	else
		loglik_star = likelihood(u_star, σ, m)
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

function sample_ν!(pars, m)

	@unpack r, ν, accept_ν, u0, a, κ, loglik, u, σ, η_0 = pars
	@unpack λ, ν_tune, ν_prior, N, T = m

	# Proposal
	forward_prop = Gamma(ν / ν_tune, ν_tune)
	ν_star = rand(forward_prop)
	back_prop = Gamma(ν_star / ν_tune, ν_tune)

	# Proposal process model
	p_star =  DEparams(r, a, κ, ν_star, λ)
	u0_star = logistic.(η_0) .* r ./ ν_star
	u_star = process_all(p_star, u0_star, m)

	# Proposal likelihood
	if size(u_star, 1) < T
		loglik_star = -Inf
	else
		loglik_star = likelihood(u_star, σ, m)
	end

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(ν_prior, ν_star) + logpdf(back_prop, ν)
	mh2 = loglik + logpdf(ν_prior, ν) + logpdf(forward_prop, ν_star)

	# Accept/reject
	prob = exp(mh1 - mh2)
	if rand() > prob
		accept_ν = 0
	else
		accept_ν = 1
		ν = ν_star
		u = u_star
		u0 = u0_star
		loglik = loglik_star
	end

	@pack! pars = ν, accept_ν, u, u0, loglik
end

# Sample functional response saturation constant
function sample_u0!(pars, m)

	@unpack r, accept_u0, u0, η_0, a, κ, loglik, u, σ, ν = pars
	@unpack λ, η_0_prior, u0_tune, N, T = m

	# Proposal
	forward_prop = MvNormal(η_0, u0_tune)
	η_0_star = rand(forward_prop)
	u0_star = logistic.(η_0_star) .* r ./ ν

	# Proposal process model
	p =  DEparams(r, a, κ, ν, λ)
	u_star = process_all(p, u0_star, m)

	# Proposal likelihood
	if size(u_star, 1) < T
		loglik_star = -Inf
	else
		loglik_star = likelihood(u_star, σ, m)
	end

	# Computing the MH ratio
	mh1 = loglik_star + logpdf(η_0_prior, η_0_star)
	mh2 = loglik + logpdf(η_0_prior, η_0)

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
	@unpack X_r, Ω_r, Ω_β_r, μ_β_r = m

	# Sample β_r
	A = Symmetric(X_r' * Ω_r * X_r + Ω_β_r)
	A_inv = inv(A)

	b = X_r' * Ω_r * η_r + Ω_β_r * μ_β_r

	β_r = rand(MvNormal(A_inv * b, A_inv))

	@pack! pars = β_r

end

# Conjugate Gibbs updates of regression coefficients on K
function sample_β_a!(pars, m)

	@unpack η_a = pars
	@unpack X_a, Ω_a, Ω_β_a, μ_β_a = m

	# Sample β_r
	A = Symmetric(X_a' * Ω_a * X_a + Ω_β_a)
	A_inv = inv(A)

	b = X_a' * Ω_a * η_a + Ω_β_a * μ_β_a

	β_a = rand(MvNormal(A_inv * b, A_inv))

	@pack! pars = β_a

end

function mcmc(m, pars, keep_every, nburn, nmcmc)

	nkeep = floor(Int64, nmcmc / keep_every)

	chain = Dict(:r => fill(0.0, m.N, nkeep),
	             :a => fill(0.0, m.N, nkeep),
				 :nu => fill(0.0, nkeep),
				 :kappa => fill(0.0, nkeep),
				 :sigma => fill(0.0, nkeep),
				 :u0 => fill(0.0, m.N, nkeep),
				 :beta_r => fill(0.0, size(m.X_r, 2), nkeep),
				 :beta_a => fill(0.0, size(m.X_a, 2), nkeep),
				 :eta_r => fill(0.0, m.N, nkeep),
				 :eta_a => fill(0.0, m.N, nkeep),
				 :eta_0 => fill(0.0, m.N, nkeep),
	             :accept_r => fill(0, nkeep),
				 :accept_a => fill(0, nkeep),
				 :accept_kappa => fill(0, nkeep),
				 :accept_nu => fill(0, nkeep),
				 :accept_sigma => fill(0, nkeep),
				 :accept_u0 => fill(0, nkeep),
				 :u => fill(0.0, m.T, m.N, nkeep),
				 :zpred => fill(0.0, m.T, m.N, nkeep))

	# Initialize process and likelihood
	p = DEparams(pars.r, pars.a, pars.κ, pars.ν, m.λ)
	pars.u = process_all(p, pars.u0, m)
	pars.loglik = likelihood(pars.u, pars.σ, m)

	# Burn-in
	for i in 1:nburn

		if i % 500 == 0
			@show i
		end

		# Sampling

		sample_ν!(pars, m)

		sample_a!(pars, m)

		sample_κ!(pars, m)

		sample_u0!(pars, m)

		sample_r!(pars, m)

		sample_σ!(pars, m)

	end

	for i in 1:nmcmc

		if i % 500 == 0
			@show i
		end

		# Sampling
		sample_ν!(pars, m)

		sample_a!(pars, m)

		sample_β_a!(pars, m)

		sample_κ!(pars, m)

		sample_u0!(pars, m)

		sample_r!(pars, m)

		sample_β_r!(pars, m)

		sample_σ!(pars, m)

		sample_z!(pars, m)

		# Saving samples
		if i % keep_every == 0
			idx = floor(Int64, i / keep_every)

			chain[:r][:, idx] = pars.r
			chain[:a][:, idx] = pars.a
			chain[:nu][idx] = pars.ν
			chain[:kappa][idx] = pars.κ
			chain[:sigma][idx] = pars.σ
			chain[:u0][:, idx] = pars.u0
			chain[:beta_r][:, idx] = pars.β_r
			chain[:beta_a][:, idx] = pars.β_a
			chain[:eta_r][:, idx] = pars.η_r
			chain[:eta_a][:, idx] = pars.η_a
			chain[:eta_0][:, idx] = pars.η_0

			chain[:accept_r][idx] = pars.accept_r
			chain[:accept_a][idx] = pars.accept_a
			chain[:accept_kappa][idx] = pars.accept_κ
			chain[:accept_nu][idx] = pars.accept_ν
			chain[:accept_sigma][idx] = pars.accept_σ
			chain[:accept_u0][idx] = pars.accept_u0

			chain[:u][:, :, idx] = pars.u
			chain[:zpred][:, :, idx] = pars.z
		end
	end

	return chain

end
