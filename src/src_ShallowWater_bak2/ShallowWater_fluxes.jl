function HyperbolicFlux!(model::SWE,
    x::Vector{MFloat}, u::Vector{MFloat}, 
    ComputeJ::Bool, f::Matrix{MFloat}, df_du::Array{MFloat,3})where MFloat<:AbstractMatrix{Float64}

    #Get variables:
    eta             = u[1]
    q1              = u[2]
    q2              = u[3]
    b               = model.b(x)
    h               = eta-b
    g               = model.g

    #Velocities:
    v1              = @tturbo @. q1/h
    v2              = @tturbo @. q2/h

    #Hyperbolic fluxes:
    @tturbo @. f[1,1]           += q1
    @tturbo @. f[1,2]           += q2
    @tturbo @. f[2,1]           += v1*q1 + g*(0.5*eta*eta-eta*b)
    @tturbo @. f[2,2]           += v2*q1
    @tturbo @. f[3,1]           += v1*q2
    @tturbo @. f[3,2]           += v2*q2 + g*(0.5*eta*eta-eta*b)

    #Get Jacobian, if necessary:
    if ComputeJ

        #Derivatives of f(1,1)
        @tturbo @. df_du[1,1,1]     += 0
        @tturbo @. df_du[1,1,2]     += 1
        @tturbo @. df_du[1,1,3]     += 0

        #Derivatives of f(1,2)
        @tturbo @. df_du[1,2,1]     += 0
        @tturbo @. df_du[1,2,2]     += 0
        @tturbo @. df_du[1,2,3]     += 1

        #Derivatives of f(2,1)
        @tturbo @. df_du[2,1,1]     += -q1^2/h^2 + g*(eta-b)
        @tturbo @. df_du[2,1,2]     += 2*q1/h
        @tturbo @. df_du[2,1,3]     += 0

        #Derivatives of f(2,2)
        @tturbo @. df_du[2,2,1]     += -q1*q2/h^2
        @tturbo @. df_du[2,2,2]     += q2/h
        @tturbo @. df_du[2,2,3]     += q1/h

        #Derivatives of f(3,1)
        @tturbo @. df_du[3,1,1]     += -q1*q2/h^2
        @tturbo @. df_du[3,1,2]     += q2/h
        @tturbo @. df_du[3,1,3]     += q1/h

        #Derivatives of f(3,2)
        @tturbo @. df_du[3,2,1]     += -q2^2/h^2 + g*(eta-b)
        @tturbo @. df_du[3,2,2]     += 0
        @tturbo @. df_du[3,2,3]     += 2*q2/h

    end

    return

end

function Source!(model::SWE, x::Vector{MFloat},
    u::Vector{MFloat}, du_dx::Matrix{MFloat}, ComputeJ::Bool,
    Q::Vector{MFloat}, dQ_du::Matrix{MFloat}, dQ_du_dx::Array{MFloat,3}) where MFloat<:AbstractMatrix{Float64}


    #Get variables:
    eta             = u[1]
    q1              = u[2]
    q2              = u[3]
    b               = model.b(x)
    db_dx, db_dy    = model.gradb(x)
    h               = eta-b
    g               = model.g
    gamma           = model.gamma

    @tturbo @. Q[1]             += 0
    @tturbo @. Q[2]             += -gamma*q1 - g*eta*db_dx
    @tturbo @. Q[3]             += -gamma*q2 - g*eta*db_dy

    if ComputeJ

        #Derivatives of Q(2)
        @tturbo @. dQ_du[2,1]       += -g*db_dx
        @tturbo @. dQ_du[2,2]       += -gamma

        #Derivatives of Q(3)
        @tturbo @. dQ_du[3,1]       += -g*db_dy
        @tturbo @. dQ_du[3,2]       += 0
        @tturbo @. dQ_du[3,3]       += -gamma
        
    end

end

function epsilonFlux!(model::SWE, tau::MFloat, duB::Matrix{MFloat},
                      ComputeJ::Bool,
                      fB::Matrix{MFloat}, dfB_dduB::Array{MFloat,4};
                      IIv::Vector{Int}=Vector{Int}(1:model.nVars)) where MFloat<:Matrix{Float64}

    for II=IIv, jj=1:2
        @avxt @. fB[II,jj]              -= tau*duB[II,jj]
    end

    if ComputeJ
        for II=IIv, jj=1:2
            @avxt @. dfB_dduB[II,jj,II,jj]  -= tau
        end
    end

    return

end
