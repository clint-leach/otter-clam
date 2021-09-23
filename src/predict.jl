function predict(chain, m)

    nmcmc = length(chain[:a])
    npred = size(m.X_r_all, 1)

    # Building storage arrays
    log_r_pred = fill(0.0, npred)
    log_a_pred = fill(0.0, npred)
    log_u_pred = fill(0.0, T, npred)
    log_flux_pred = fill(0.0, T, npred)

    Σ_pred_r = PDMat(Symmetric(m.Σuu_r - m.Σuo_r * m.Ω_r * m.Σuo_r'))
    Σ_pred_a = PDMat(Symmetric(m.Σuu_a - m.Σuo_a * m.Ω_a * m.Σuo_a'))
    Σ_pred_0 = PDMat(Symmetric(m.Σuu_0 - m.Σuo_0 * m.Ω_0 * m.Σuo_0'))

    for i in 1:1

        if i % 500 == 0
			@show i
		end

        β_r = chain[:beta_r][:, i]
        β_a = chain[:beta_a][:, i]
        η_r = chain[:eta_r][:, i]
        η_a = chain[:eta_a][:, i]
        η_0 = chain[:eta_0][:, i]
        κ = chain[:kappa][i]
        ν = chain[:nu][i]

        # Compute μ_pred for r and K
        μ_pred_r = m.X_r_all * β_r + m.Σuo_r * m.Ω_r * (η_r - m.X_r * β_r)
        μ_pred_a = m.X_a_all * β_a + m.Σuo_a * m.Ω_a * (η_a - m.X_a * β_a)
        μ_pred_0 = fill(m.μ_η_0, npred) + m.Σuo_0 * m.Ω_0 * (η_0 - fill(m.μ_η_0, m.N))

        # Draw r and K from predictive distributions
        η_r_pred = rand(MvNormal(μ_pred_r, Σ_pred_r))
        η_a_pred = rand(MvNormal(μ_pred_a, Σ_pred_a))
        η_0_pred = rand(MvNormal(μ_pred_0, Σ_pred_0))

        r = exp.(η_r_pred)
        a = exp.(η_a_pred)
        u0 = logistic.(η_0_pred) .* r ./ ν
        
        # Run process model at every site
        p =  DEparams(r, a, κ, ν, λ_all)
        u = process_pred(p, u0, m)

        flux = similar(u)
        for j in 1:npred
            for t in 1:T
                flux[t, j] = a[j] * u[t, j] ^ 2 * λ_all[j](t)  / (u[t, j] ^ 2 + κ ^ 2)
            end
        end

        # Compute running means of r, u, and flux
        log_r_pred = log_r_pred + 1.0 / nmcmc * η_r_pred
        log_a_pred = log_a_pred + 1.0 / nmcmc * η_a_pred
        log_u_pred = log_u_pred + 1.0 / nmcmc * log.(u)
        log_flux_pred = log_flux_pred + 1.0 / nmcmc * log.(flux)
    end

    preds = Dict(:log_r => log_r_pred,
                 :log_a => log_a_pred,
                 :log_u => log_u_pred,
                 :log_flux => log_flux_pred)

    return preds
end


