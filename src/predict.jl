function predict(chain, iters, m)

    @unpack X_all, X, Σuo_r, Ω_r, Σuu_r, Σuo_0, Ω_0, Σuu_0, λ_all, T = m

    # nmcmc = length(chain["sigma"])
    nmcmc = length(iters)
    npred = size(X_all, 1)

    # Building storage arrays
    r_pred = fill(0.0, npred, nmcmc)
    u_pred = fill(0.0, T, npred, nmcmc)

    Σ_pred_r = PDMat(Symmetric(Σuu_r - Σuo_r * Ω_r * Σuo_r'))
    Σ_pred_0 = PDMat(Symmetric(Σuu_0 - Σuo_0 * Ω_0 * Σuo_0'))

    for i in 1:nmcmc

        @show i

        β_r = chain["beta_r"][:, iters[i]]
        β_0 = chain["beta_0"][:, iters[i]]
        η_r = chain["eta_r"][:, iters[i]]
        η_0 = chain["eta_0"][:, iters[i]]
        a = chain["a"][iters[i]]
        κ = chain["kappa"][iters[i]]
        K = chain["K"][iters[i]]

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

        # Return r and K vectors and U matrix
        r_pred[:, i] = r
        u_pred[:, :, i] = u

    end

    preds = Dict("r" => r_pred,
                 "u" => u_pred)

    return preds
end


