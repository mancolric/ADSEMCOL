function HyperbolicFlux!(model::NHWW,
    u::Vector{MFloat},
    ComputeJ::Bool, f::Matrix{MFloat}, df_du::Array{MFloat,3}) where MFloat<:AbstractMatrix{Float64}

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
    g               = model.g
    c               = model.c
    eta             = u[1]
    q1              = u[2]
    q2              = u[3]
    q3              = u[4]
    P               = u[5]
    b               = u[6]
    h               = eta-b
    xi_h            = @tturbo @. P/(h*h)
    p               = @tturbo @. (c*c)/3.0 * xi_h * (1.0-xi_h)
    
    #Hyperbolic fluxes:
    @tturbo @. f[1,1]           += q1
    @tturbo @. f[1,2]           += q2
    @tturbo @. f[2,1]           += q1*q1/h + g*(0.5*eta*eta-eta*b) + h*p
    @tturbo @. f[2,2]           += q1*q2/h
    @tturbo @. f[3,1]           += q2*q1/h
    @tturbo @. f[3,2]           += q2*q2/h + g*(0.5*eta*eta-eta*b) + h*p
    @tturbo @. f[4,1]           += q3*q1/h
    @tturbo @. f[4,2]           += q3*q2/h
    @tturbo @. f[5,1]           += P*q1/h
    @tturbo @. f[5,2]           += P*q2/h

    #Get Jacobian, if necessary. Note that
    #   h       = eta - b
    #   xi_h    = P/h^2
    #   p       = 1/3 c^2 (xi_h - (xi_h)^2)
    #Therefore, terms containing "h" have to be derived w.r.t. eta and b, and terms
    #containing "p" have to be derived w.r.t. h (and therefore eta and b) and P.
    #In particular:
    #   dh/deta     = 1
    #   dh/db       = -1
    #   dp/deta     = dp/d(xi_h) * dxi_h/dh * 1     = 1/3 c^2 (1 - 2 xi_h) * (-2P)/h^3
    #   dp/db       = dp/d(xi_h) * dxi_h/dh * (-1)  = 1/3 c^2 (1 - 2 xi_h) * (+2P)/h^3
    #   dp/dP       = dp/d(xi_h) * dxi_h/dP         = 1/3 c^2 (1 - 2 xi_h) * 1/h^2
    if ComputeJ

        #Auxiliary derivatives:
        dh_deta                     = 1.0
        dh_db                       = -1.0
        dp_deta                     = @tturbo @. (c*c)/3*(1-2*xi_h)*(-2*P)/(h*h*h)
        dp_db                       = @tturbo @. -dp_deta
        dp_dP                       = @tturbo @. (c*c)/3*(1-2*xi_h)/(h*h)
        
        #Derivatives of f(1,1)
        @tturbo @. df_du[1,1,2]     += 1.0

        #Derivatives of f(1,2)
        @tturbo @. df_du[1,2,3]     += 1.0

        #Derivatives of f(2,1)
        @tturbo @. df_du[2,1,1]     += -q1^2/h^2*dh_deta + g*(eta-b) + dh_deta*p + h*dp_deta
        @tturbo @. df_du[2,1,2]     += 2*q1/h
        @tturbo @. df_du[2,1,5]     += h*dp_dP
        @tturbo @. df_du[2,1,6]     += -q1^2/h^2*dh_db - g*eta + dh_db*p + h*dp_db
        
        #Derivatives of f(2,2)
        @tturbo @. df_du[2,2,1]     += -q1*q2/h^2*dh_deta
        @tturbo @. df_du[2,2,2]     += q2/h
        @tturbo @. df_du[2,2,3]     += q1/h
        @tturbo @. df_du[2,2,6]     += -q1*q2/h^2*dh_db

        #Derivatives of f(3,1)
        @tturbo @. df_du[3,1,1]     += -q1*q2/h^2*dh_deta
        @tturbo @. df_du[3,1,2]     += q2/h
        @tturbo @. df_du[3,1,3]     += q1/h
        @tturbo @. df_du[3,1,6]     += -q1*q2/h^2*dh_db
        
        #Derivatives of f(3,2)
        @tturbo @. df_du[3,2,1]     += -q2^2/h^2*dh_deta + g*(eta-b) + dh_deta*p + h*dp_deta
        @tturbo @. df_du[3,2,3]     += 2*q2/h
        @tturbo @. df_du[3,2,5]     += h*dp_dP
        @tturbo @. df_du[3,2,6]     += -q2^2/h^2*dh_db - g*eta + dh_db*p + h*dp_db

        #Derivatives of f(4,1)
        @tturbo @. df_du[4,1,1]     += -q3*q1/h^2*dh_deta
        @tturbo @. df_du[4,1,2]     += q3/h
        @tturbo @. df_du[4,1,4]     += q1/h
        @tturbo @. df_du[4,1,6]     += -q3*q1/h^2*dh_db
        
        #Derivatives of f(4,2)
        @tturbo @. df_du[4,2,1]     += -q3*q2/h^2*dh_deta
        @tturbo @. df_du[4,2,3]     += q3/h
        @tturbo @. df_du[4,2,4]     += q2/h
        @tturbo @. df_du[4,2,6]     += -q3*q2/h^2*dh_db
        
        #Derivatives of f(5,1)
        @tturbo @. df_du[5,1,1]     += -P*q1/h^2*dh_deta
        @tturbo @. df_du[5,1,2]     += P/h
        @tturbo @. df_du[5,1,5]     += q1/h
        @tturbo @. df_du[5,1,6]     += -P*q1/h^2*dh_db
        
        #Derivatives of f(5,2)
        @tturbo @. df_du[5,2,1]     += -P*q2/h^2*dh_deta
        @tturbo @. df_du[5,2,3]     += P/h
        @tturbo @. df_du[5,2,5]     += q2/h
        @tturbo @. df_du[5,2,6]     += -P*q2/h^2*dh_db

    end

    return

end

function Source1!(model::NHWW, x::Vector{MFloat}, tau_char::Float64,
    u::Vector{MFloat}, du_dx::Matrix{MFloat}, ComputeJ::Bool,
    Q::Vector{MFloat}, dQ_du::Matrix{MFloat}, dQ_du_dx::Array{MFloat,3}) where MFloat<:AbstractMatrix{Float64}


    #Get variables:
    g               = model.g
    c               = model.c
    eta             = u[1]
    q1              = u[2]
    q2              = u[3]
    q3              = u[4]
    P               = u[5]
    b               = u[6]
    h               = eta-b
    xi_h            = @tturbo @. P/(h*h)
    p               = @tturbo @. (c*c)/3.0 * xi_h * (1.0-xi_h)
    b_exact         = model.b(x)
    db_dx           = du_dx[6,1]
    db_dy           = du_dx[6,2]

    #Source terms:
    @tturbo @. Q[2]             += - g*eta*db_dx - 3/2*p/xi_h*db_dx
    @tturbo @. Q[3]             += - g*eta*db_dy - 3/2*p/xi_h*db_dy
    @tturbo @. Q[4]             += (c*c)*(1.0-xi_h)
    @tturbo @. Q[5]             += q3 - 3/2*q1*db_dx - 3/2*q2*db_dy
    @tturbo @. Q[6]             += -1.0/tau_char*(b-b_exact)

    #Get Jacobian, if necessary. Note that
    #   h       = eta - b
    #   xi_h    = P/h^2
    #   p       = 1/3 c^2 (xi_h - (xi_h)^2)
    #Therefore, terms containing "h" have to be derived w.r.t. eta and b,
    #and terms containing "xi_h" or "p" have to be derived w.r.t. eta, b and P.
    #In particular:
    #   dh/deta     = 1
    #   dh/db       = -1
    #   dxi_h/deta  = dxi_h/dh * 1      = -2P/h^3
    #   dxi_h/db    = dxi_h/dh * (-1)   = 2P/h^3
    #   dxi_h/dP    = 1/h^2
    #   dp/deta     = dp/d(xi_h) * dxi_h/dh * 1     = 1/3 c^2 (1 - 2 xi_h) * (-2P)/h^3
    #   dp/db       = dp/d(xi_h) * dxi_h/dh * (-1)  = 1/3 c^2 (1 - 2 xi_h) * (+2P)/h^3
    #   dp/dP       = dp/d(xi_h) * dxi_h/dP         = 1/3 c^2 (1 - 2 xi_h) * 1/h^2
    if ComputeJ

        #Auxiliary derivatives:
        dh_deta                     = 1.0
        dh_db                       = -1.0
        dxih_deta                   = @tturbo @. -2*P/(h*h*h)
        dxih_db                     = @tturbo @. -dxih_deta
        dxih_dP                     = @tturbo @. 1.0/(h*h)
        dp_deta                     = @tturbo @. (c*c)/3*(1-2*xi_h)*(-2*P)/(h*h*h)
        dp_db                       = @tturbo @. -dp_deta
        dp_dP                       = @tturbo @. (c*c)/3*(1-2*xi_h)/(h*h)
        
        #Derivatives of Q(2)
        @tturbo @. dQ_du[2,1]       += -g*db_dx - 3/2*dp_deta/xi_h*db_dx + 3/2*p/xi_h^2*dxih_deta*db_dx
        @tturbo @. dQ_du[2,5]       += - 3/2*dp_dP/xi_h*db_dx + 3/2*p/xi_h^2*dxih_dP*db_dx
        @tturbo @. dQ_du[2,6]       += - 3/2*dp_db/xi_h*db_dx + 3/2*p/xi_h^2*dxih_db*db_dx
        @tturbo @. dQ_du_dx[2,6,1]  += -g*eta - 3/2*p/xi_h

        #Derivatives of Q(3)
        @tturbo @. dQ_du[3,1]       += -g*db_dy - 3/2*dp_deta/xi_h*db_dy + 3/2*p/xi_h^2*dxih_deta*db_dy
        @tturbo @. dQ_du[3,5]       += - 3/2*dp_dP/xi_h*db_dy + 3/2*p/xi_h^2*dxih_dP*db_dy
        @tturbo @. dQ_du[3,6]       += - 3/2*dp_db/xi_h*db_dy + 3/2*p/xi_h^2*dxih_db*db_dy
        @tturbo @. dQ_du_dx[3,6,2]  += -g*eta - 3/2*p/xi_h
        
        #Derivatives of Q(4)
        @tturbo @. dQ_du[4,1]       += -(c*c)*dxih_deta
        @tturbo @. dQ_du[4,5]       += -(c*c)*dxih_dP
        @tturbo @. dQ_du[4,6]       += -(c*c)*dxih_db
        
        #Derivatives of Q(5)
        @tturbo @. dQ_du[5,2]       += -3/2*db_dx
        @tturbo @. dQ_du[5,3]       += -3/2*db_dy
        @tturbo @. dQ_du[5,4]       += 1.0
        @tturbo @. dQ_du_dx[5,6,1]  += -3/2*q1
        @tturbo @. dQ_du_dx[5,6,2]  += -3/2*q2
        
        #Derivatives of Q(6)
        @tturbo @. dQ_du[6,6]       += -1.0/tau_char
        
    end
    
    return

end
