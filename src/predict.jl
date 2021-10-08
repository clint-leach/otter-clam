function predict(chain, m)

    nmcmc = length(chain[:nu])
    npred = size(m.X_r_all, 1)

    # Building storage arrays
    r_pred = fill(0.0, npred, nmcmc)
    a_pred = fill(0.0, npred, nmcmc)
    u_pred = fill(0.0, T, npred, nmcmc)
    flux_pred = fill(0.0, T, npred, nmcmc)

    Σ_pred_r = PDMat(Symmetric(m.Σuu_r - m.Σuo_r * m.Ω_r * m.Σuo_r'))
    Σ_pred_a = PDMat(Symmetric(m.Σuu_a - m.Σuo_a * m.Ω_a * m.Σuo_a'))

    for i in 1:nmcmc

        if i % 500 == 0
			@show i
		end

        β_r = chain[:beta_r][:, i]
        β_a = chain[:beta_a][:, i]
        η_r = chain[:eta_r][:, i]
        η_a = chain[:eta_a][:, i]
        κ = chain[:kappa][i]
        ν = chain[:nu][i]

        # Compute μ_pred for r and K
        μ_pred_r = m.X_r_all * β_r + m.Σuo_r * m.Ω_r * (η_r - m.X_r * β_r)
        μ_pred_a = m.X_a_all * β_a + m.Σuo_a * m.Ω_a * (η_a - m.X_a * β_a)

        # Draw r and K from predictive distributions
        η_r_pred = rand(MvNormal(μ_pred_r, Σ_pred_r))
        η_a_pred = rand(MvNormal(μ_pred_a, Σ_pred_a))

        r = [η_r_pred[j] > 0 ? η_r_pred[j] : 0.0 for j in 1:npred]
        a = [η_a_pred[j] > 0 ? η_a_pred[j] : 0.0 for j in 1:npred]
        u0 = logistic(4.0) .* r ./ ν
        
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
        r_pred[:, i] = r #r_pred + 1.0 / nmcmc * r
        a_pred[:, i] = a #a_pred + 1.0 / nmcmc * a
        u_pred[:, :, i] = u #log_u_pred + 1.0 / nmcmc * log.(u)
        flux_pred[:, :, i] = flux #flux_pred + 1.0 / nmcmc * flux
    end

    preds = Dict(:r => r_pred,
                 :a => a_pred,
                 :u => u_pred,
                 :flux => flux_pred)

    return preds
end


