@warn "Review calc_cpi, calc_cvi, calc_h"

#-----------------------------------------------------------------------------------------------------
# Function to compute dp_du:
function calc_dp_du(model::GasH2, u::Vector{Matrix{Float64}}, T::Matrix{Float64}, dT_drhoYi::Vector{Matrix{Float64}}, dT_dmi::Vector{Matrix{Float64}}, dT_dE::Matrix{Float64})

    R           = model.R
    Wi          = model.Wi
    nSpecies    = model.nSpecies
    nVars       = model.nVars
    dp_du       = Vector{Matrix{Float64}}(undef,nVars)

    prod_rhoY   = zeros(size(u[1]))

    for j = 1:nSpecies
        @. prod_rhoY += u[j]/Wi[j]
    end

    # dp_drhoYi
    for i = 1:nSpecies
       dp_du[i]     = @. R*T/Wi[i] + dT_drhoYi[i]*prod_rhoY*R
    end
    # dp_dmi
    dp_du[nVars-2]  = dT_dmi[1].*prod_rhoY*R
    dp_du[nVars-1]  = dT_dmi[2].*prod_rhoY*R
    # dp_dE
    dp_du[nVars]    = dT_dE.*prod_rhoY*R

    return dp_du
end

#----------------------------------------------------------------------------------------

# evaluate cp and cv from NASA polynomials coefficients:
function pol_c(coef::Vector{Float64}, T::VecOrMat{Float64})

    a = reverse(coef)
    p = a[1].*ones(size(T))
    for i = 2:length(a)
        @. p = p*T + a[i]
    end

    return p
end

# evaluate h from NASA polynomials coefficients
function pol_h(coef::Vector{Float64}, T::VecOrMat{Float64})

    a = reverse(coef)
    b = length(a)
    p = a[1]/b.*ones(size(T))
    for i = 2:length(a)
        b = b - 1
        @. p = p*T + a[i]/b
    end

    return p
end

#-----------------------------------------------------------------------------------------------------
# cpi
function calc_cpi(model::GasH2, T::Matrix{Float64})
    N       = model.nSpecies
    cpi     = Vector{Matrix{Float64}}(undef, N)
    for ii=1:N
        cpi[ii]   = zeros(size(T))
    end
    coef1   = model.coef1
    coef2   = model.coef2
    coef3   = model.coef3
    Wi      = model.Wi

    # cpi/R = a1/T^2 + a2/T + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
    for i = 1:N
        # 200-1000 K
       aux     = @. 200<=T<1000
       cpi[i][aux]  = @. (model.R/Wi[i])*(coef1[i,1]/T[aux]^2 + coef1[i,2]/T[aux] + $pol_c(coef1[i,3:7],T[aux]))
#        cpi[i][aux]  = @. (model.R/Wi[i])*(coef1[i,1]/T[aux]^2 + coef1[i,2]/T[aux] + coef1[i,3] + coef1[i,4]*T[aux] + coef1[i,5]*T[aux]^2 + coef1[i,6]*T[aux]^3 + coef1[i,7]*T[aux]^4)

       # 1000-6000 K
       aux     = @. 1000<=T<6000
       cpi[i][aux]  = @. (model.R/Wi[i])*(coef2[i,1]/T[aux]^2 + coef2[i,2]/T[aux] + $pol_c(coef2[i,3:7],T[aux]))
#        cpi[i][aux]  = @. (model.R/Wi[i])*(coef2[i,1]/T[aux]^2 + coef2[i,2]/T[aux] + coef2[i,3] + coef2[i,4]*T[aux] + coef2[i,5]*T[aux]^2 + coef2[i,6]*T[aux]^3 + coef2[i,7]*T[aux]^4)

       # 6000-20000 K
       aux     = @. 6000<=T<=20000
       cpi[i][aux]  = @. (model.R/Wi[i])*(coef3[i,1]/T[aux]^2 + coef3[i,2]/T[aux] + $pol_c(coef3[i,3:7],T[aux]))
#        cpi[i][aux]  = @. (model.R/Wi[i])*(coef3[i,1]/T[aux]^2 + coef3[i,2]/T[aux] + coef3[i,3] + coef3[i,4]*T[aux] + coef3[i,5]*T[aux]^2 + coef3[i,6]*T[aux]^3 + coef3[i,7]*T[aux]^4)
    end

    return cpi
end
# cp
function calc_cp_total(model::GasH2, T::Matrix{Float64}, Yi::Vector{Matrix{Float64}})
    N       = model.nSpecies
    cp_total = zeros(size(T))
    cpi      = calc_cpi(model,T)

    for i = 1:N
        @. cp_total += Yi[i]*cpi[i]
    end

    return cp_total
end
#-----------------------------------------------------------------------------------------------------
# cvi
function calc_cvi(model::GasH2, T::Matrix{Float64})
    N       = model.nSpecies
    cvi     = Vector{Matrix{Float64}}(undef, N)
    for ii=1:N
        cvi[ii]   = zeros(size(T))
    end
    coef1   = model.coef1
    coef2   = model.coef2
    coef3   = model.coef3
    Wi      = model.Wi

    # cpi/R = a1/T^2 + a2/T + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
    for i = 1:N
        # 200-1000 K
        aux     = @. 200<=T<1000
        cvi[i][aux]  = @. (model.R/Wi[i])*(coef1[i,1]/T[aux]^2 + coef1[i,2]/T[aux] + $pol_c(coef1[i,3:7],T[aux]) - 1.0)
#         cvi[i][aux]  = @. (model.R/Wi[i])*(coef1[i,1]/T[aux]^2 + coef1[i,2]/T[aux] + coef1[i,3] + coef1[i,4]*T[aux] + coef1[i,5]*T[aux]^2 + coef1[i,6]*T[aux]^3 + coef1[i,7]*T[aux]^4 - 1.0)

        # 1000-6000 K
        aux     = @. 1000<=T<6000
        cvi[i][aux]  = @. (model.R/Wi[i])*(coef2[i,1]/T[aux]^2 + coef2[i,2]/T[aux] + $pol_c(coef2[i,3:7],T[aux]) - 1.0)
#         cvi[i][aux]  = @. (model.R/Wi[i])*(coef2[i,1]/T[aux]^2 + coef2[i,2]/T[aux] + coef2[i,3] + coef2[i,4]*T[aux] + coef2[i,5]*T[aux]^2 + coef2[i,6]*T[aux]^3 + coef2[i,7]*T[aux]^4 - 1.0)

        # 6000-20000 K
        aux     = @. 6000<=T<=20000
        cvi[i][aux]  = @. (model.R/Wi[i])*(coef3[i,1]/T[aux]^2 + coef3[i,2]/T[aux] + $pol_c(coef3[i,3:7],T[aux]) - 1.0)
#         cvi[i][aux]  = @. (model.R/Wi[i])*(coef3[i,1]/T[aux]^2 + coef3[i,2]/T[aux] + coef3[i,3] + coef3[i,4]*T[aux] + coef3[i,5]*T[aux]^2 + coef3[i,6]*T[aux]^3 + coef3[i,7]*T[aux]^4 - 1.0)
    end

    return cvi
end
# cv
function calc_cv_total(model::GasH2, T::Matrix{Float64}, Yi::Vector{Matrix{Float64}})
    N       = model.nSpecies
    cv_total = zeros(size(T))
    cvi      = calc_cvi(model,T)

    for i = 1:N
        @. cv_total += Yi[i]*cvi[i]
    end

    return cv_total
end
#-----------------------------------------------------------------------------------------------------
# hi
function calc_hi(model::GasH2, T::Matrix{Float64})
    N       = model.nSpecies
    hi      = Vector{Matrix{Float64}}(undef, N)
    for ii=1:N
        hi[ii]   = zeros(size(T))
    end
    coef1   = model.coef1
    coef2   = model.coef2
    coef3   = model.coef3
    Wi      = model.Wi

    for i = 1:N
        # 200-1000 K
        aux     = @. 200<=T<1000
        hi[i][aux]    = @. model.R/Wi[i]*(-coef1[i,1]/T[aux] + coef1[i,2]*log(T[aux]) + $pol_h(coef1[i,3:7],T[aux])*T[aux] + coef1[i,8])
#         hi[i][aux]    = @. model.R/Wi[i]*(-coef1[i,1]/T[aux] + coef1[i,2]*log(T[aux]) + coef1[i,3]*T[aux] + (coef1[i,4]/2)*T[aux]^2 + (coef1[i,5]/3)*T[aux]^3 + (coef1[i,6]/4)*T[aux]^4 + (coef1[i,7]/5)*T[aux]^5 + coef1[i,8])

        # 1000-6000 K
        aux     = @. 1000<=T<6000
        hi[i][aux]    = @. model.R/Wi[i]*(-coef2[i,1]/T[aux] + coef2[i,2]*log(T[aux]) + $pol_h(coef2[i,3:7],T[aux])*T[aux] + coef2[i,8])
#         hi[i][aux]    = @. model.R/Wi[i]*(-coef2[i,1]/T[aux] + coef2[i,2]*log(T[aux]) + coef2[i,3]*T[aux] + (coef2[i,4]/2)*T[aux]^2 + (coef2[i,5]/3)*T[aux]^3 + (coef2[i,6]/4)*T[aux]^4 + (coef2[i,7]/5)*T[aux]^5 + coef2[i,8])

        # 6000-20000 K
        aux     = @. 6000<=T<=20000
        hi[i][aux]    = @. model.R/Wi[i]*(-coef3[i,1]/T[aux] + coef3[i,2]*log(T[aux]) + $pol_h(coef3[i,3:7],T[aux])*T[aux] + coef3[i,8])
#         hi[i][aux]    = @. model.R/Wi[i]*(-coef3[i,1]/T[aux] + coef3[i,2]*log(T[aux]) + coef3[i,3]*T[aux] + (coef3[i,4]/2)*T[aux]^2 + (coef3[i,5]/3)*T[aux]^3 + (coef3[i,6]/4)*T[aux]^4 + (coef3[i,7]/5)*T[aux]^5 + coef3[i,8])
    end

    return hi
end
# h
function calc_h_total(model::GasH2, T::Matrix{Float64}, Yi::Vector{Matrix{Float64}})
    N       = model.nSpecies
    h_total = zeros(size(T))
    hi      = calc_hi(model,T)

    for i = 1:N
        @. h_total += Yi[i]*hi[i]
    end

    return h_total
end
#-----------------------------------------------------------------------------------------------------
# ei
function calc_ei(model::GasH2, T::Matrix{Float64})
    N  = model.nSpecies
    ei = Vector{Matrix{Float64}}(undef, N)
    Wi = model.Wi
    hi = calc_hi(model, T)

    for i = 1:N
        ei[i] = @. hi[i] - (model.R/Wi[i])*(T-model.T0)
    end

    return ei
end
# e
function calc_e_total(model::GasH2, T::Matrix{Float64}, Yi::Vector{Matrix{Float64}})
    N       = model.nSpecies
    e_total = zeros(size(T))
    ei      = calc_ei(model,T)

    for i = 1:N
        @. e_total += Yi[i]*ei[i]
    end

    return e_total
end
#-----------------------------------------------------------------------------------------------------
# Find T with Newton-Raphson method:
function calc_T(model::GasH2, emat::Matrix{Float64}, Yi::Vector{Matrix{Float64}})

    T_k     = fill(model.T0, size(e))
    e_T     = zeros(size(e))
    cv_T    = zeros(size(e))
    T_pk    = zeros(size(e))

    # Newton-Raphson
    for n = 1:10

        # inicializamos e y cv
        @. e_T      = 0.0
        @. cv_T     = 0.0

        # e(Tk)
        e_T     = calc_e_total(model,T_k,Yi)
        # cv(Tk)
        cv_T    = calc_cv_total(model,T_k,Yi)

        # Algoritmo
        @. T_pk = abs(T_k - (e_T - emat)/cv_T)

        # tolerancia
        if norm(T_pk - T_k, Inf) < 1e-8
            break
        elseif n==10
            @warn "Reached maximum of iterations when computing T"
        end

        # update Tk
        @. T_k = T_pk
    end

#         figure()
#         scatter(T, T.-T_pk)
#         title("Err T-Tpk")

    return T_pk
    
end

#-----------------------------------------------------------------------------
#DEBUG:

# comprobación de resultados con diferencias finitas
function test_fun_aux_H2()

    var = crear_T_Yi()
    T   = var[1]
    Yi  = var[2]
## Comprobación e y c

    # e(T)
    e_T     = calc_e_total(GasH2(),T,Yi)
    figure()
    plot(T[:],e_T[:], "b")
    title("e(T)")

    # cv(T)
    cv_1  = calc_cv_total(GasH2(),T,Yi)
    figure()
    plot(T[:],cv_1[:], "b")

    # ddff cv(T)
    delta   = 1e-3
    cvi_2   = Vector{Matrix{Float64}}(undef,length(Yi))
    for ii=1:length(Yi)
        cvi_2[ii]   = zeros(size(T))
    end
    cv_2  = zeros(size(T))

    for j = 1:length(T)
        T_pert      = deepcopy(T)
        T_pert[1,j] = T[1,j] - delta
        e_neg       = calc_ei(GasH2(),T_pert)
        T_pert[1,j] = T[1,j] + delta
        e_pos       = calc_ei(GasH2(),T_pert)
        for i = 1:length(Yi)
            cvi_2[i][j] = (e_pos[i][1,j]-e_neg[i][1,j])/(2*delta)
        end
    end

    for i = 1:length(Yi)
        @. cv_2 += Yi[i]*cvi_2[i]
    end

    plot(T,cv_2, "r")
    title("cv(T)")

    figure()
    scatter(T,cv_1.-cv_2)
    title("Err cv1-cv2")

## Comprobación h y cp

    # h(T)
    h_T     = calc_h_total(GasH2(),T,Yi)
    figure()
    plot(T[:],h_T[:])
    title("h(T)")

    # cp(T)
    cp_1  = calc_cp_total(GasH2(),T,Yi)
    figure()
    plot(T[:],cp_1[:])

    # ddff cp(T)
    delta   = 1e-3
    cpi_2   = Vector{Matrix{Float64}}(undef,length(Yi))
    for ii=1:length(Yi)
        cpi_2[ii]   = zeros(size(T))
    end
    cp_2  = zeros(size(T))

    for j = 1:length(T)
        T_pert      = deepcopy(T)
        T_pert[1,j] = T[1,j] - delta
        h_neg       = calc_hi(GasH2(),T_pert)
        T_pert[1,j] = T[1,j] + delta
        h_pos       = calc_hi(GasH2(),T_pert)
        for i = 1:length(Yi)
            cpi_2[i][j] = (h_pos[i][1,j]-h_neg[i][1,j])/(2*delta)
        end

    end

    for i = 1:length(Yi)
        @. cp_2 += Yi[i]*cpi_2[i]
    end

    plot(T,cp_2)
    title("cp(T)")

    figure()
    scatter(T,cp_1.-cp_2)
    title("Err cp1-cp2")

end

# funcion para generar T y Yi y evaluar
function crear_T_Yi()

    N       = 1000
    T       = zeros(1,N); T[1:end] = LinRange(300, 800, N)
    Yi      = [fill(rand(),1,N), fill(rand(),1,N), fill(rand(),1,N), fill(rand(),1,N), fill(rand(),1,N),
               fill(rand(),1,N), fill(rand(),1,N), fill(rand(),1,N), fill(rand(),1,N), fill(rand(),1,N), fill(rand(),1,N),fill(rand(),1,N), fill(rand(),1,N)]

    return T,Yi
end
