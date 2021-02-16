# Model object storing all constants

@with_kw struct model

	z
	T
	tobs
	tspan
	N
	λ

	r_tune
	r_prior

	a_tune
	a_prior

	κ_tune
	κ_prior

	K_tune
	K_prior

	u0_tune
	u0_prior
end

@with_kw mutable struct parameters

	u0::Vector{Float64}
	r::Float64
	K::Float64
	a::Float64
	κ::Float64

	u
	loglik::Float64

	accept_r::Int64
	accept_a::Int64
	accept_κ::Int64
	accept_K::Int64
	accept_u0::Int64

end
