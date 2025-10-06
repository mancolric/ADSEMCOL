function HyperbolicFlux!(model::NHWW,
    u::Vector{MFloat},
    ComputeJ::Bool, f::Matrix{MFloat}, df_du::Array{MFloat,3})where MFloat<:AbstractMatrix{Float64}

    return HyperbolicFlux1!(model, u, ComputeJ, f, df_du)

end

function Source!(model::NHWW, x::Vector{MFloat}, tau_char::Float64,
    u::Vector{MFloat}, du_dx::Matrix{MFloat}, ComputeJ::Bool,
    Q::Vector{MFloat}, dQ_du::Matrix{MFloat}, dQ_du_dx::Array{MFloat,3}) where MFloat<:AbstractMatrix{Float64}

    return Source1!(model, x, tau_char, u, du_dx, ComputeJ, Q, dQ_du, dQ_du_dx)

end

function epsilonFlux!(model::NHWW, tau::MFloat, duB::Matrix{MFloat},
                      ComputeJ::Bool,
                      fB::Matrix{MFloat}, dfB_dduB::Array{MFloat,4};
                      IIv::Vector{Int}=Vector{Int}(1:model.nVars-1)) where MFloat<:Matrix{Float64}

    for II=IIv, jj=1:2
        @avxt @. fB[II,jj]              -= tau*duB[II,jj]
    end

    if ComputeJ
        for II=IIv, jj=1:2
            @avxt @. dfB_dduB[II,jj,II,jj]      -= tau
        end
    end

    return

end

#------------------------------

function HyperbolicFlux1!(model::NHWW,
    u::Vector{MFloat},
    ComputeJ::Bool, f::Matrix{MFloat}, df_du::Array{MFloat,3})where MFloat<:AbstractMatrix{Float64}

    #Get variables:
    eta             = u[1]
    q1              = u[2]
    q2              = u[3]
    q3              = u[4]
    p               = u[5]
    b               = u[6]
    h               = eta-b
    g               = model.g
    gamma           = model.gamma
    c               = model.c
    
    #Hyperbolic fluxes:
    @tturbo @. f[1,1]           += q1
    @tturbo @. f[1,2]           += q2
    @tturbo @. f[2,1]           += q1*q1/h + g*(0.5*eta*eta-eta*b) + h*p
    @tturbo @. f[2,2]           += q1*q2/h
    @tturbo @. f[3,1]           += q2*q1/h
    @tturbo @. f[3,2]           += q2*q2/h + g*(0.5*eta*eta-eta*b) + h*p
    @tturbo @. f[4,1]           += q3*q1/h
    @tturbo @. f[4,2]           += q3*q2/h
    @tturbo @. f[5,1]           += p*q1/h + c*c*q1/h
    @tturbo @. f[5,2]           += p*q2/h + c*c*q2/h

    #Get Jacobian, if necessary. Note that d/deta = d/dh * 1, d/db = d/dh * (-1), since h=eta-b
    if ComputeJ

        #Derivatives of f(1,1)
        @tturbo @. df_du[1,1,2]     += 1.0

        #Derivatives of f(1,2)
        @tturbo @. df_du[1,2,3]     += 1.0

        #Derivatives of f(2,1)
        @tturbo @. df_du[2,1,1]     += -q1^2/h^2 + g*(eta-b) + p
        @tturbo @. df_du[2,1,2]     += 2*q1/h
        @tturbo @. df_du[2,1,5]     += h
        @tturbo @. df_du[2,1,6]     += +q1^2/h^2 - g*eta - p
        
        #Derivatives of f(2,2)
        @tturbo @. df_du[2,2,1]     += -q1*q2/h^2
        @tturbo @. df_du[2,2,2]     += q2/h
        @tturbo @. df_du[2,2,3]     += q1/h
        @tturbo @. df_du[2,2,6]     += +q1*q2/h^2

        #Derivatives of f(3,1)
        @tturbo @. df_du[3,1,1]     += -q1*q2/h^2
        @tturbo @. df_du[3,1,2]     += q2/h
        @tturbo @. df_du[3,1,3]     += q1/h
        @tturbo @. df_du[3,1,6]     += +q1*q2/h^2
        
        #Derivatives of f(3,2)
        @tturbo @. df_du[3,2,1]     += -q2^2/h^2 + g*(eta-b) + p
        @tturbo @. df_du[3,2,3]     += 2*q2/h
        @tturbo @. df_du[3,2,5]     += h
        @tturbo @. df_du[3,2,6]     += +q2^2/h^2 - g*eta - p

        #Derivatives of f(4,1)
        @tturbo @. df_du[4,1,1]     += -q3*q1/h^2
        @tturbo @. df_du[4,1,2]     += q3/h
        @tturbo @. df_du[4,1,4]     += q1/h
        @tturbo @. df_du[4,1,6]     += +q3*q1/h^2
        
        #Derivatives of f(4,2)
        @tturbo @. df_du[4,2,1]     += -q3*q2/h^2
        @tturbo @. df_du[4,2,3]     += q3/h
        @tturbo @. df_du[4,2,4]     += q2/h
        @tturbo @. df_du[4,2,6]     += +q3*q2/h^2
        
        #Derivatives of f(5,1)
        @tturbo @. df_du[5,1,1]     += -p*q1/h^2 - c*c*q1/h^2
        @tturbo @. df_du[5,1,2]     += p/h + c*c/h
        @tturbo @. df_du[5,1,5]     += q1/h
        @tturbo @. df_du[5,1,6]     += +p*q1/h^2 + c*c*q1/h^2
        
        #Derivatives of f(5,2)
        @tturbo @. df_du[5,2,1]     += -p*q2/h^2 - c*c*q2/h^2
        @tturbo @. df_du[5,2,3]     += p/h + c*c/h
        @tturbo @. df_du[5,2,5]     += q2/h
        @tturbo @. df_du[5,2,6]     += +p*q2/h^2 + c*c*q2/h^2

    end

    return

end

function Source1!(model::NHWW, x::Vector{MFloat}, tau_char::Float64,
    u::Vector{MFloat}, du_dx::Matrix{MFloat}, ComputeJ::Bool,
    Q::Vector{MFloat}, dQ_du::Matrix{MFloat}, dQ_du_dx::Array{MFloat,3}) where MFloat<:AbstractMatrix{Float64}


    #Get variables:
    eta             = u[1]
    q1              = u[2]
    q2              = u[3]
    q3              = u[4]
    p               = u[5]
    b               = u[6]
    h               = eta-b
    b_exact         = model.b(x)
    g               = model.g
    gamma           = model.gamma
    c               = model.c
    db_dx           = du_dx[6,1]
    db_dy           = du_dx[6,2]

    @tturbo @. Q[2]             += - g*eta*db_dx - gamma*p*db_dx
    @tturbo @. Q[3]             += - g*eta*db_dy - gamma*p*db_dy
    @tturbo @. Q[4]             += gamma*p
    @tturbo @. Q[5]             += -(p+c*c)*(2*q3/h^2 - 2*q1/h^2*db_dx - 2*q2/h^2*db_dy) 
    @tturbo @. Q[6]             += -1.0/tau_char*(b-b_exact)

    if ComputeJ

        #Derivatives of Q(2)
        @tturbo @. dQ_du[2,1]       += -g*db_dx
        @tturbo @. dQ_du[2,5]       += -gamma*db_dx
        @tturbo @. dQ_du_dx[2,6,1]  += -g*eta - gamma*p

        #Derivatives of Q(3)
        @tturbo @. dQ_du[3,1]       += -g*db_dy
        @tturbo @. dQ_du[3,5]       += -gamma*db_dy
        @tturbo @. dQ_du_dx[3,6,2]  += -g*eta - gamma*p
        
        #Derivatives of Q(4)
        @tturbo @. dQ_du[4,5]       += gamma
        
        #Derivatives of Q(5)
        @tturbo @. dQ_du[5,1]       += -(p+c*c)*(-4*q3/h^3 + 4*q1/h^3*db_dx + 4*q2/h^3*db_dy)
        @tturbo @. dQ_du[5,2]       += -(p+c*c)*(-2/h^2*db_dx)
        @tturbo @. dQ_du[5,3]       += -(p+c*c)*(-2/h^2*db_dy)
        @tturbo @. dQ_du[5,4]       += -(p+c*c)*(2/h^2)
        @tturbo @. dQ_du[5,5]       += -(2*q3/h^2 - 2*q1/h^2*db_dx - 2*q2/h^2*db_dy) 
        @tturbo @. dQ_du[5,6]       += +(p+c*c)*(-4*q3/h^3 + 4*q1/h^3*db_dx + 4*q2/h^3*db_dy)
        @tturbo @. dQ_du_dx[5,6,1]  += +(p+c*c)*(2*q1/h^2)
        @tturbo @. dQ_du_dx[5,6,2]  += +(p+c*c)*(2*q2/h^2)
        
        #Derivatives of Q(6)
        @tturbo @. dQ_du[6,6]       += -1.0/tau_char
        
    end
    
    return

end
