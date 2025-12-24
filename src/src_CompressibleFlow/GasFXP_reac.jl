function calc_mdot(model::GasFXP, rhoY::Vector{Matrix{Float64}}, RT::Matrix{Float64})

    #Arrhenius coeffs:
    kI          = @tturbo @. model.BI*exp(-model.RTI/RT)
    kB          = @tturbo @. model.BB*exp(-model.RTB/RT)
    kR          = fill(model.BR, size(RT))
    
    #mdot:
    mdot        = Vector{Matrix{Float64}}(undef, model.nSpecies)
    mdot[1]     = @tturbo @. - rhoY[1]*kI - rhoY[1]*rhoY[2]*kB
    mdot[2]     = @tturbo @.   rhoY[1]*kI + rhoY[1]*rhoY[2]*kB - rhoY[2]*kR
    mdot[3]     = @tturbo @.                                   + rhoY[2]*kR
    
    return mdot
    
end

function calc_dmdot(model::GasFXP, rhoY::Vector{Matrix{Float64}}, RT::Matrix{Float64})
    
    delta       = 1e-6;
    N           = model.nSpecies
    dfun_drhoY  = Matrix{Matrix{Float64}}(undef, N, N)
    dfun_dRT    = Vector{Matrix{Float64}}(undef, N)

    # dfun/drhoY
    for j = 1:N
    
        #Copy vector of pointers. Deep copy rhoY[j]:
        rhoY_pert                   = copy(rhoY)
        rhoY_pert[j]                = copy(rhoY[j])
        
        #Perturb rhoY and evaluate:
        @tturbo @. rhoY_pert[j]     = rhoY[j]-delta
        f_neg                       = calc_mdot(model, rhoY_pert, RT)
        @tturbo @. rhoY_pert[j]     = rhoY[j]+delta
        f_pos                       = calc_mdot(model, rhoY_pert, RT)
        
        # Finite differences:
        for i=1:N
            dfun_drhoY[i,j] = (f_pos[i]-f_neg[i])/(2*delta)
        end
    
    end

    # dfun/dT

    RT_pert     = @tturbo @. RT-delta
    f_neg       = calc_mdot(model, rhoY, RT_pert)
    RT_pert     = @tturbo @. RT+delta
    f_pos       = calc_mdot(model, rhoY, RT_pert)
    # Finite differences:
    for i=1:N
        dfun_dRT[i]     = (f_pos[i]-f_neg[i])/(2*delta)
    end

    return dfun_drhoY, dfun_dRT
    
end
