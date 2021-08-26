function predict(chain, m)

    @unpack X_all, X, Σuo_r, Ω_r, Σuu_r, Σuo_0, Ω_0, Σuu_0, λ_all, T = m

    nmcmc = length(chain["a"])
    npred = size(X_all, 1)

    # Building storage arrays
    log_r_pred = fill(0.0, npred)
    log_u_pred = fill(0.0, T, npred)
    log_flux_pred = fill(0.0, T, npred)

    Σ_pred_r = PDMat(Symmetric(Σuu_r - Σuo_r * Ω_r * Σuo_r'))
    Σ_pred_0 = PDMat(Symmetric(Σuu_0 - Σuo_0 * Ω_0 * Σuo_0'))

    for i in 1:nmcmc

        if i % 500 == 0
			@show i
		end

        β_r = chain[:beta_r][:, i]
        β_0 = chain[:beta_0][:, i]
        η_r = chain[:eta_r][:, i]
        η_0 = chain[:eta_0][:, i]
        a = chain[:a][i]
        κ = chain[:kappa][i]
        K = chain[:K][i]

        # Compute μ_pred for r and K
        μ_pred_r = X_all * β_r + Σuo_r * Ω_r * (η_r - X * β_r)
        μ_pred_0 = X_all * β_0 + Σuo_0 * Ω_0 * (η_0 - X * β_0)

        # Draw r and K from predictive distributions
        η_r_pred = rand(MvNormal(μ_pred_r, Σ_pred_r))
        η_0_pred = rand(MvNormal(μ_pred_0, Σ_pred_0))

        r = exp.(η_r_pred)
        u0 = K * logistic.(η_0_pred)
        
        # Run process model at every site
        p =  DEparams(r, a, κ, K, λ_all)
        u = process_all(p, u0, m)

        flux = similar(u)
        for j in 1:npred
            for t in 1:T
                u[t, j] = u[t, j] > 0 ? u[t, j] : 1e-10
                flux[t, j] = a * u[t, j] * λ_all[j](t)  / (u[t, j] + κ)
            end
        end

        # Compute running means of r, u, and flux
        log_r_pred = log_r_pred + 1.0 / nmcmc * η_r_pred
        log_u_pred = log_u_pred + 1.0 / nmcmc * log.(u)
        log_flux_pred = log_flux_pred + 1.0 / nmcmc * log.(flux)
    end

    preds = Dict(:log_r => log_r_pred,
                 :log_u => log_u_pred,
                 :log_flux => log_flux_pred)

    return preds
end


