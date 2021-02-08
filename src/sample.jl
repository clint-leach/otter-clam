# Sampling scripts and helper functions

function likelihood(u, m)

	loglik = 0.0
	for t in 1:m.T
		loglik += logpdf(truncated(Normal(u[t], 5.0), 0.0, Inf), z[t])
	end

	return loglik
end

# Sample prey population growth rate
function sample_r!(pars, m)

	@unpack r, accept_r, u0, K, a, κ, loglik, u = pars
	@unpack λ, r_tune, r_prior = m

	r_tune = 0.0001

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

function mcmc(m, pars, nmcmc)

	chain = Dict("r" => fill(0.0, nmcmc),
	             "accept_r" => fill(0, nmcmc),
				 "u" => fill(0.0, 26, nmcmc))

	# Initialize process and likelihood
	p =  [pars.r, pars.K, pars.a, pars.κ, m.λ]
	pars.u = process(p, pars.u0, m)
	pars.loglik = likelihood(pars.u, m)

	@progress for i in 1:nmcmc

		sample_r!(pars, m)

		chain["r"][i] = pars.r
		chain["accept_r"][i] = pars.accept_r
		chain["u"][:, i] = pars.u

	end

	return chain

end
