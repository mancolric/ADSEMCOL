function calc_mdot(model::GasFXP, rhoY::Vector{Matrix{Float64}}, RT::Matrix{Float64})

    #Arrhenius coeffs:
    kI          = @tturbo @. model.BI*exp(-model.RTI/RT)
    kB          = @tturbo @. model.BB*exp(-model.RTB/RT)
    kR          = fill(model.BR, size(RT))
    
    #mdot:
    mdot        = Vector{Matrix{Float64}}(undef, model.nSpecies)
    mdot[1]     = - rhoY[1]*kI - rhoY[1]*rhoY[2]*kB
    mdot[2]     =   rhoY[1]*kI + rhoY[1]*rhoY[2]*kB - rhoY[3]*kR
    mdot[3]     =                                   + rhoY[3]*kR
    
    return mdot
    
end

function calc_dmdot(model::GasFXP, rhoY::Vector{Matrix{Float64}}, T::Matrix{Float64})
    
    delta       = 1e-6;
    N           = model.nSpecies
    dfun_drhoY  = Matrix{Matrix{Float64}}(undef, N, N)
    dfun_dT     = Vector{Matrix{Float64}}(undef, N)

    # dfun/drhoY
    for j = 1:N
    
        #Copy vector of pointers. Deep copy rhoY[j]:
        rhoY_pert                   = copy(rhoY)
        rhoY_pert[j]                = copy(rhoY[j])
        
        #Perturb rhoY and evaluate:
        @tturbo @. rhoY_pert[j]     = rhoY[j]-delta
        f_neg                       = calc_mdot(model, rhoY_pert, T)
        @tturbo @. rhoY_pert[j]     = rhoY[j]+delta
        f_pos                       = calc_mdot(model, rhoY_pert, T)
        
        # Finite differences:
        for i=1:N
            dfun_drhoY[i,j] = (f_pos[i]-f_neg[i])/(2*delta)
        end
    
    end

    # dfun/dT

    T_pert      = @tturbo @. T-delta
    f_neg       = calc_mdot(model, rhoY, T_pert)
    T_pert      = @tturbo @. T+delta
    f_pos       = calc_mdot(model, rhoY, T_pert)
    # Finite differences:
    for i=1:N
        dfun_dT[i] = (f_pos[i]-f_neg[i])/(2*delta)
    end

    return dfun_drhoY, dfun_dT
    
end
