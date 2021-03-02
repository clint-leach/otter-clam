# Model object storing all constants

@with_kw struct model

	z
	T
	tobs
	tspan
	N
	λ

	α_r_tune
	α_r_prior
	μ_r
	L_r

	a_tune
	a_prior

	κ_tune
	κ_prior

	α_K_tune
	α_K_prior
	μ_K
	L_K

	u0_tune
	u0_prior

	σ_tune
	σ_prior
end

@with_kw mutable struct parameters

	u0::Vector{Float64}
	log_r::Vector{Float64}
	α_r::Vector{Float64}
	log_K::Vector{Float64}
	α_K::Vector{Float64}
	a::Float64
	κ::Float64
	σ::Float64

	u::Matrix{Float64}
	loglik::Vector{Float64}

	accept_r::Int64
	accept_a::Int64
	accept_κ::Int64
	accept_K::Int64
	accept_u0::Vector{Int64}
	accept_σ::Int64

end

@with_kw struct DEparams

	log_r::Vector{Float64}
	log_K::Vector{Float64}
	a::Float64
	κ::Float64
	λ::Array{Interpolations.BSplineInterpolation{Float64,1,Array{Float64,1},BSpline{Linear},Tuple{Base.OneTo{Int64}}},1}

end
