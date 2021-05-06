function predict(chain, m)

    @unpack X_all, X, Σuo_r, Ω_r, Σuu_r, Σuo_K, Ω_K, Σuu_K, λ_all, T = m

    nmcmc = length(chain["sigma"])
    npred = size(X_all, 1)

    # Building storage arrays
    r_pred = fill(0.0, npred, nmcmc)
    K_pred = fill(0.0, npred, nmcmc)
    u_pred = fill(0.0, T, npred, nmcmc)

    Σ_pred_r = PDMat(Symmetric(Σuu_r - Σuo_r * Ω_r * Σuo_r'))
    Σ_pred_K = PDMat(Symmetric(Σuu_K - Σuo_K * Ω_K * Σuo_K'))

    for i in 1:nmcmc

        @show i

        β_r = chain["beta_r"][:, i]
        β_K = chain["beta_K"][:, i]
        log_r = chain["r"][:, i]
        log_K = chain["K"][:, i]
        a = chain["a"][i]
        κ = chain["kappa"][i]

        # Compute μ_pred for r and K
        μ_pred_r = X_all * β_r + Σuo_r * Ω_r * (log_r - X * β_r)
        μ_pred_K = X_all * β_K + Σuo_K * Ω_K * (log_K - X * β_K)

        # Draw r and K from predictive distributions
        log_r_pred = rand(MvNormal(μ_pred_r, Σ_pred_r))
        log_K_pred = rand(MvNormal(μ_pred_K, Σ_pred_K))

        # Run process model at every site
        u0 = log_K_pred
        p =  DEparams(log_r_pred, log_K_pred, a, κ, λ_all)
        u = process_all(p, u0, m)

        # Return r and K vectors and U matrix
        r_pred[:, i] = log_r_pred
        K_pred[:, i] = log_K_pred
        u_pred[:, :, i] = u

    end

    preds = Dict("r" => r_pred,
                 "K" => K_pred,
                 "u" => u_pred)

    return preds
end


