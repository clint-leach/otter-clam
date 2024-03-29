# Model object storing all constants

@with_kw struct model

	# Data and constants
	z
	λ
	X_r
	X_a
	quads
	nq = size(z, 1)
	T = size(z, 2)
	N = size(z, 3)
	tspan = (1.0, Float64(T))

	# Scalar priors
	ν_tune
	ν_prior

	σ_tune
	σ_prior

	# Distance matrices
	Doo

	# r regression priors
	Ω_β_r
	μ_β_r

	# r spatial priors
	σ_r
	ρ_r
	Σ_r = PDMat(σ_r ^ 2 * exp.(-0.5 * (Doo ./ ρ_r)))
	Ω_r = inv(Σ_r)
	A_r = Symmetric(Ω_r - Ω_r * X_r * inv(X_r' * Ω_r * X_r + Ω_β_r) * X_r' * Ω_r)
	b_r = μ_β_r' * Ω_β_r * inv(X_r' * Ω_r * X_r + Ω_β_r) * X_r' * Ω_r
	η_r_prior = MvNormal(inv(A_r) * b_r', inv(A_r))

	r_tune

	# a regression priors
	Ω_β_a
	μ_β_a

	# a spatial priors
	σ_a
	ρ_a
	Σ_a = PDMat(σ_a ^ 2 * exp.(-0.5 * (Doo ./ ρ_a)))
	Ω_a = inv(Σ_a)
	A_a = Symmetric(Ω_a - Ω_a * X_a * inv(X_a' * Ω_a * X_a + Ω_β_a) * X_a' * Ω_a)
	b_a = μ_β_a' * Ω_β_a * inv(X_a' * Ω_a * X_a + Ω_β_a) * X_a' * Ω_a
	η_a_prior = MvNormal(inv(A_a) * b_a', inv(A_a))

	a_tune
	
end

@with_kw mutable struct parameters

	r::Vector{Float64}
	η_r::Vector{Float64} = log.(r)
	η_r_pred::Vector{Float64} = η_r
	β_r::Vector{Float64}
	β_end::Vector{Float64} = β_r

	a::Vector{Float64}
	η_a::Vector{Float64} = log.(a)
	η_a_pred::Vector{Float64} = η_a
	β_a::Vector{Float64}

	ν::Float64
	σ::Vector{Float64}

	u0::Vector{Float64} = r
	u::Matrix{Float64}
	z::Matrix{Float64}
	loglik::Float64 = 0.0
	z_mean::Matrix{Float64} = similar(z)
	z_var::Matrix{Float64} = similar(z)

	accept_r::Int64 = 0
	accept_a::Int64 = 0
	accept_ν::Int64 = 0
	accept_σ::Int64 = 0

end

@with_kw struct DEparams

	r::Vector{Float64}
	a::Vector{Float64}
	ν::Float64
	λ::Vector{Interpolations.BSplineInterpolation{Float64, 1, Vector{Float64}, BSpline{Linear{Throw{OnGrid}}}, Tuple{Base.OneTo{Int64}}}} 

end
