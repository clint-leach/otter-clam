# Model object storing all constants

@with_kw struct model

	z
	T
	tspan
	λ

	r_tune
	r_prior

end

@with_kw mutable struct parameters

	u0::Float64
	r::Float64
	K::Float64
	a::Float64
	κ::Float64

	u::Vector{Float64}
	loglik::Float64

	accept_r::Int64

end
