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
	K_tune
	K_prior

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
	η_r_prior = MvNormal(inv(A_r) * b_r', inv(A_r))

	r_tune

	# u0 regression priors
	Ω_β_0
	μ_β_0

	# K spatial priors
	σ_0
	ρ_0
	Σ_0 = PDMat(σ_0 ^ 2 * exp.(-0.5 * (Doo ./ ρ_0)))
	Ω_0 = inv(Σ_0)
	A_0 = Symmetric(Ω_0 - Ω_0 * X * inv(X' * Ω_0 * X + Ω_β_0) * X' * Ω_0)
	b_0 = μ_β_0' * Ω_β_0 * inv(X' * Ω_0 * X + Ω_β_0) * X' * Ω_0
	η_0_prior = MvNormal(inv(A_0) * b_0', inv(A_0))

	u0_tune
	
	# Prediction components
	X_all
	λ_all
	Σuo_r = σ_r ^ 2 * exp.(-0.5 * (Duo ./ ρ_r))
	Σuu_r = PDMat(σ_r ^ 2 * exp.(-0.5 * (Duu ./ ρ_r)))
	Σuo_0 = σ_0 ^ 2 * exp.(-0.5 * (Duo ./ ρ_0))
	Σuu_0 = PDMat(σ_0 ^ 2 * exp.(-0.5 * (Duu ./ ρ_0)))

end

@with_kw mutable struct parameters

	u0::Vector{Float64}
	η_0::Vector{Float64}
	β_0::Vector{Float64}

	r::Vector{Float64}
	η_r::Vector{Float64}
	β_r::Vector{Float64}

	K::Float64
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
	accept_u0::Int64
	accept_σ::Int64

end

@with_kw struct DEparams

	r::Vector{Float64}
	a::Float64
	κ::Float64
	K::Float64
	λ::Array{Interpolations.BSplineInterpolation{Float64,1,Array{Float64,1},BSpline{Linear},Tuple{Base.OneTo{Int64}}},1}

end
