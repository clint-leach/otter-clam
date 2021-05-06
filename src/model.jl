# Model object storing all constants

@with_kw struct model

	# Data and constants
	z
	λ
	X
	T = size(z, 1)
	N = size(z, 2)
	tspan = (1.0, Float64(T))
	p = size(X, 2)

	# Scalar priors
	a_tune
	a_prior

	κ_tune
	κ_prior

	σ_tune
	σ_prior

	# Distance matrices
	Doo
	Duo
	Duu

	# r regression priors
	Ω_β_r
	μ_β_r

	# r spatial priors
	σ_r
	ρ_r
	Σ_r = PDMat(σ_r ^ 2 * exp.(-0.5 * (Doo ./ ρ_r)))
	Ω_r = inv(Σ_r)
	A_r = Symmetric(Ω_r - Ω_r * X * inv(X' * Ω_r * X + Ω_β_r) * X' * Ω_r)
	b_r = μ_β_r' * Ω_β_r * inv(X' * Ω_r * X + Ω_β_r) * X' * Ω_r
	log_r_prior = MvNormal(inv(A_r) * b_r', inv(A_r))

	r_tune

	# K regression priors
	Ω_β_K
	μ_β_K

	# K spatial priors
	σ_K
	ρ_K
	Σ_K = PDMat(σ_K ^ 2 * exp.(-0.5 * (Doo ./ ρ_K)))
	Ω_K = inv(Σ_K)
	A_K = Symmetric(Ω_K - Ω_K * X * inv(X' * Ω_K * X + Ω_β_K) * X' * Ω_K)
	b_K = μ_β_K' * Ω_β_K * inv(X' * Ω_K * X + Ω_β_K) * X' * Ω_K
	log_K_prior = MvNormal(inv(A_K) * b_K', inv(A_K))

	K_tune
	
	# Prediction components
	X_all
	λ_all
	Σuo_r = σ_r ^ 2 * exp.(-0.5 * (Duo ./ ρ_r))
	Σuu_r = PDMat(σ_r ^ 2 * exp.(-0.5 * (Duu ./ ρ_r)))
	Σuo_K = σ_K ^ 2 * exp.(-0.5 * (Duo ./ ρ_K))
	Σuu_K = PDMat(σ_K ^ 2 * exp.(-0.5 * (Duu ./ ρ_K)))

end

@with_kw mutable struct parameters

	u0::Vector{Float64}

	log_r::Vector{Float64}
	β_r::Vector{Float64}

	log_K::Vector{Float64}
	β_K::Vector{Float64}

	a::Float64
	κ::Float64
	σ::Float64

	u::Matrix{Float64}
	loglik::Vector{Float64}

	accept_r::Int64
	accept_a::Int64
	accept_κ::Int64
	accept_K::Int64
	accept_β::Int64
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
