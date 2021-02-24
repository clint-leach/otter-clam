# Model object storing all constants

@with_kw struct model

	z
	T
	tobs
	tspan
	N
	λ

	Σ_r_tune
	log_r_prior

	a_tune
	a_prior

	κ_tune
	κ_prior

	Σ_K_tune
	log_K_prior

	u0_tune
	u0_prior

	σ_tune
	σ_prior
end

@with_kw mutable struct parameters

	u0::Vector{Float64}
	log_r::Vector{Float64}
	log_K::Vector{Float64}
	a::Float64
	κ::Float64
	σ::Float64

	u
	loglik::Float64

	accept_r::Int64
	accept_a::Int64
	accept_κ::Int64
	accept_K::Int64
	accept_u0::Int64
	accept_σ::Int64

end
