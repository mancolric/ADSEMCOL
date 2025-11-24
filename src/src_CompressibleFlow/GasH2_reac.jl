# Arrhenius and falloff K functions:
function K_arrhenius(B::Float64, n::Float64, Ta::Float64, T::Matrix{Float64})

    K = @. B*exp(n*log(T)-Ta/T)

    return K
end

function K_falloff(B::Vector{Float64}, n::Vector{Float64}, Ta::Vector{Float64}, T::Matrix{Float64}, Fc::Float64, Cm::Matrix{Float64})

    K0   = K_arrhenius(B[1],n[1],Ta[1],T)
    Kinf = K_arrhenius(B[2],n[2],Ta[2],T)
    Pr   = @. K0*Cm/Kinf

    C = -0.4-0.67*log10(Fc)
    N = 0.75-1.27*log10(Fc)

    f = @. (1 + C/log10(Pr))/((N-0.14*C)/log10(Pr) - 0.14) # Expresión más robusta
    F = @. Fc^(1/(1+f^2))

    K = @. (Kinf*Pr*F)/(1+Pr)

    return K
end

#-------------------------------------------------------------------------------
# computation of mdot
function calc_mdot(model::GasH2, u::Vector{Matrix{Float64}}, T::Matrix{Float64})

    N   = model.nSpecies
    Wi  = model.Wi
    C   = Vector{Matrix{Float64}}(undef, N)
    for ii=1:N
        C[ii] = zeros(size(u[1]))
    end
    # N species
    for i = 1:N
        @. C[i] = 1e-6*u[i]/Wi[i] # kg/m^3 / kg/mol -> mol/cm^3
    end

    ## Concentrations------------------------------------------
    # Species
    C_H2    = C[1]
    C_O2    = C[2]
    C_H2O   = C[3]
    C_N2    = C[4]
    C_He    = C[5]
    C_Ar    = C[6]
    C_CO    = C[7]
    C_CO2   = C[8]
    C_H     = C[9]
    C_OH    = C[10]
    C_HO2   = C[11]
    C_H2O2  = C[12]
    C_O     = C[13]

    # Chaperon coeffs
    A       = model.coefChaperon
    # Indefinite molecules
    C_Ma    = sum(C.*A[1,:])
    C_Mb    = sum(C.*A[2,:])
    C_Mc    = sum(C.*A[3,:])
    C_Md    = sum(C.*A[4,:])

#     C_Ma    = @. 2.5*C_H2 + C_O2 + 16.0*C_H2O + C_N2 + 0.7*(C_He + C_Ar) + 1.2*C_CO + 2.4*C_CO2 + C_H + C_OH + C_HO2 + C_H2O2 + C_O
#     C_Mb    = @. 2.5*C_H2 + C_O2 + 12.0*C_H2O + C_N2 + 0.4*(C_He + C_Ar) + 1.9*C_CO + 3.8*C_CO2 + C_H + C_OH + C_HO2 + C_H2O2 + C_O
#     C_Mc    = @. 2.5*C_H2 + C_O2 + 12.0*C_H2O + C_N2 + 0.5*(C_He + C_Ar) + 1.9*C_CO + 3.8*C_CO2 + C_H + C_OH + C_HO2 + C_H2O2 + C_O
#     C_Md    = @. 2.5*C_H2 + C_O2 + 6.0*C_H2O + C_N2 + 0.4*C_He + 0.7*C_Ar + 1.5*C_CO + 2.0*C_CO2 + C_H + C_OH + C_HO2 + 6.0*C_H2O2 + C_O

    ## Kj ------------------------------------------------------
    Q   = model.nReactions
    K   = Vector{Matrix{Float64}}(undef,Q)
    C_Q = Vector{Vector{Matrix{Float64}}}(undef,Q)

#Q1
#1f # H + O2 → OH + O
    B           = 3.52e16
    n           = -0.7
    Ta          = 8590.0
    K[1]        = K_arrhenius(B, n, Ta, T)

    C_Q[1]  = [C_H, C_O2]

#Q2
#1b # OH + O → H + O2
    B           = 1.05e14
    n           = -0.313
    Ta          = 132.0
    K[2]        = K_arrhenius(B, n, Ta, T)

    C_Q[2]  = [C_OH, C_O]

#Q3
#2f # H2 + O → OH + H
    B           = 5.06e4
    n           = 2.67
    Ta          = 3165.0
    K[3]        = K_arrhenius(B, n, Ta, T)

    C_Q[3]  = [C_H2, C_O]

#Q4
#2b # OH + H → H2 + O
    B           = 2.94e4
    n           = 2.64
    Ta          = 2430.0
    K[4]        = K_arrhenius(B, n, Ta, T)

    C_Q[4]  = [C_OH, C_H]

#Q5
#3f # H2 + OH → H2O +H
    B           = 1.17e9
    n           = 1.3
    Ta          = 1825.0
    K[5]        = K_arrhenius(B, n, Ta, T)

    C_Q[5]  = [C_H2, C_OH]

#Q6
#3b # H2O + H → H2 + OH
    B           = 1.42e10
    n           = 1.18
    Ta          = 9379.0
    K[6]        = K_arrhenius(B, n, Ta, T)

    C_Q[6]  = [C_H2O, C_H]

#Q7
#4f # H + O2 + Ma → HO2 + Ma
    B           = [5.75e19, 4.65e12]
    n           = [-1.4, 0.44]
    Ta          = [0.0, 0.0]
    Fc          = 0.5
    Cm          = abs.(C_Ma) # OJO posible error de puntero

    K[7]        = K_falloff(B, n, Ta, T, Fc, Cm)

    C_Q[7]  = [C_H, C_O2, C_Ma]

#Q8
#5f # HO2 + H → OH + OH
    B           = 7.08e13
    n           = 0.0
    Ta          = 148.0
    K[8]        = K_arrhenius(B, n, Ta, T)

    C_Q[8]  = [C_HO2, C_H]

#Q9
#6f # HO2 + H → H2 + O2
    B           = 1.66e13
    n           = 0.0
    Ta          = 414.0
    K[9]        = K_arrhenius(B, n, Ta, T)

    C_Q[9]  = [C_HO2, C_H]

#Q10
#6b # H2 + O2 → HO2 + H
    B           = 2.93e12
    n           = 0.356
    Ta          = 27929.0
    K[10]       = K_arrhenius(B, n, Ta, T)

    C_Q[10] = [C_H2, C_O2]

#Q11
#7f # HO2 + OH → H2O + O2
    B           = 2.89e13
    n           = 0.0
    Ta          = -250.0

    B2          = 4.5e14
    n2          = 0.0
    Ta2         = 5500.0
    K[11]       = K_arrhenius(B, n, Ta, T) + K_arrhenius(B2, n2, Ta2, T)

    C_Q[11] = [C_HO2, C_OH]

#Q12
#8f # H + OH + Mb → H2O + Mb
    B           = 4.0e22
    n           = -2.0
    Ta          = 0.0
    K[12]       = K_arrhenius(B, n, Ta, T)

    C_Q[12] = [C_H, C_OH, C_Mb]

#Q13
#8b # H2O + Mb → H + OH + Mb
    B           = 2.02e23
    n           = -1.85
    Ta          = 59682.0
    K[13]       = K_arrhenius(B, n, Ta, T)

    C_Q[13] = [C_H2O, C_Mb]

#Q14
#9f # H + H + Mc → H2 + Mc
    B           = 1.3e18
    n           = -1.0
    Ta          = 0.0
    K[14]       = K_arrhenius(B, n, Ta, T)

    C_Q[14] = [C_H, C_H, C_Mc]

#Q15
#9b # H2 + Mc → H + H + Mc
    B           = 5.41e17
    n           = -0.727
    Ta          = 52129.0
    K[15]       = K_arrhenius(B, n, Ta, T)

    C_Q[15] = [C_H2, C_Mc]

#Q16
#10f# H2O2 + Md → OH + OH + Md
    B           = [7.6e30, 2.63e19]
    n           = [-4.2, -1.27]
    Ta          = [25703.0, 25703.0]
    Fc          = 0.43
    Cm          = abs.(C_Md) # OJO posible error de puntero
    K[16]       = K_falloff(B, n, Ta, T, Fc, Cm)

    C_Q[16] = [C_H2O2, C_Md]

#Q17
#11f# HO2 + HO2 → H2O2 + O2
    B           = 1.03e14
    n           = 0.0
    Ta          = 5556.0

    B2          = 1.94e11
    n2          = 0.0
    Ta2         = -709.0
    K[17]       = K_arrhenius(B, n, Ta, T) + K_arrhenius(B2, n2, Ta2, T)

    C_Q[17] = [C_HO2, C_HO2]

#Q18
#12f# HO2 + H2 → H2O2 + H
    B           = 7.8e10
    n           = 0.61
    Ta          = 12045.0
    K[18]       = K_arrhenius(B, n, Ta, T)

    C_Q[18] = [C_HO2, C_H2]

    ## wj ------------------------------------------------------
    w       = Vector{Matrix{Float64}}(undef,Q)
    for ii=1:Q
        w[ii] = zeros(size(T))
    end

    for j=1:Q
       C_prod = 1.0
       for jj=1:length(C_Q[j])
           C_prod = C_prod.*C_Q[j][jj]
       end
       @. w[j] = K[j]*C_prod
    end

    #w is in mol/(cm^3 s). Write w in mol/(m^3 s)
    w   = w.*1e6

    ## mdot ------------------------------------------------------
    mdot       = Vector{Matrix{Float64}}(undef,N)
    for ii=1:N
        mdot[ii] = zeros(size(T))
    end

    # H2
    @. mdot[1]  = Wi[1]*(-w[3]+w[4]-w[5]+w[6]+w[9]-w[10]+w[14]-w[15]-w[18])

    # O2
    @. mdot[2]  = Wi[2]*(-w[1]+w[2]-w[7]+w[9]-w[10]+w[11]+w[17])

    # H2O
    @. mdot[3]  = Wi[3]*(w[5]-w[6]+w[11]+w[12]-w[13])

    # N2, He, Ar, CO, CO2
    @. mdot[4]= 0.0*w[1]
    @. mdot[5]= 0.0*w[1]
    @. mdot[6]= 0.0*w[1]
    @. mdot[7]= 0.0*w[1]
    @. mdot[8]= 0.0*w[1]

    # H
    @. mdot[9]  = Wi[9]*(-w[1]+w[2]+w[3]-w[4]+w[5]-w[6]-w[7]-w[8]-w[9]+w[10]-w[12]+w[13]-2*w[14]+2*w[15]+w[18])

    # OH
    @. mdot[10] = Wi[10]*(w[1]-w[2]+w[3]-w[4]-w[5]+w[6]+2*w[8]-w[11]-w[12]+w[13]+2*w[16])

    # HO2
    @. mdot[11] = Wi[11]*(w[7]-w[8]-w[9]+w[10]-w[11]-2*w[17]-w[18])

    # H2O2
    @. mdot[12] = Wi[12]*(-w[16]+w[17]+w[18])

    # O
    @. mdot[13] = Wi[13]*(w[1]-w[2]-w[3]+w[4])

    return mdot
end

#-------------------------------------------------------------------------------
# computation of dmdot
function calc_dmdot(model::GasH2, rhoY::Vector{Matrix{Float64}}, T::Matrix{Float64})
    delta       = 1e-6;
    N           = model.nSpecies
    dfun_drhoY  = Matrix{Matrix{Float64}}(undef, N, N)
    dfun_dT     = Vector{Matrix{Float64}}(undef, N)

    # dfun/drhoY
    for j = 1:N
        rhoY_pert      = deepcopy(rhoY)
        rhoY_pert[j]   = rhoY[j].-delta
        f_neg       = calc_mdot(model, rhoY_pert, T)
        rhoY_pert[j]   = rhoY[j].+delta
        f_pos       = calc_mdot(model, rhoY_pert, T)
        # Finite differences:
        for i=1:N
            dfun_drhoY[i,j] = (f_pos[i]-f_neg[i])/(2*delta)
        end
    end

    # dfun/dT

    T_pert      = deepcopy(T)
    T_pert      = T.-delta
    f_neg       = calc_mdot(model, rhoY, T_pert)
    T_pert      = T.+delta
    f_pos       = calc_mdot(model, rhoY, T_pert)
    # Finite differences:
    for i=1:N
        dfun_dT[i] = (f_pos[i]-f_neg[i])/(2*delta)
    end

    return dfun_drhoY, dfun_dT
end

#-----------------------------------------------------------------------------------------------------
# Functions for product of vectors of matrices
function prod_dmdot(dm_dT::Vector{Matrix{Float64}}, dT::Matrix{Matrix{Float64}})

    N       = length(dm_dT)
    prod    = Matrix{Matrix{Float64}}(undef,N,N)
    for i = 1:N
       for j = 1:N
           prod[i,j] = dm_dT[i].*dT[j]
       end
    end

    return prod
end

function prod_dmdot(dm_dT::Vector{Matrix{Float64}}, dT::Matrix{Float64})

    N       = length(dm_dT)
    prod    = Vector{Matrix{Float64}}(undef,N)
    for i = 1:N
        prod[i] = dm_dT[i].*dT
    end

    return prod
end

#-----------------------------------------------------------------------------
#DEBUG:

# funcion para generar u y evaluar
function test_mdot()

    n = 1
    m = 1
    u = [fill(1.2172e-01/9,n,m) for i = 1:13]
    for i = 5:8
#         u[i] = 0.0*u[i]
    end

    T = fill(1500.0,n,m);

    return calc_mdot(GasH2(),u,T)
    
end
