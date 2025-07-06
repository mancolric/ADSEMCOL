#-------------------------------------------------------------------------------
#INITIAL CONDITIONS:

u0_Gauss(x,K,sigma2)    = @. @views exp(-K*(x[:,1]^2+x[:,2]^2)/2/sigma2)

function u0_SlottedCylinder(xm,x0,y0,r0,d,w)

    #Distance to the center of the cylinder:
    rv      = @views @. sqrt((xm[:,1]-x0)^2+(xm[:,2]-y0)^2)
    uv      = zeros(size(xm,1))

    #Find points inside the cylinder:
    T       = @views @. (rv<r0) & ( (xm[:,1]>x0-r0+d) |  (abs(xm[:,2])>w/2) )
    uv[T]   .= 1

    return  uv

end

function u0_RampX(xm,x1,y1,x2,y2)

    uv      = zeros(size(xm,1))
    for ii=1:size(xm,1)
        if xm[ii,1]<=x1
            uv[ii]  = y1
        elseif xm[ii,1]>=x2
            uv[ii]  = y2
        else
            uv[ii]  = y1 + (y2-y1)/(x2-x1)*(xm[ii,1]-x1)
        end
    end

    return uv

end

function u0_hump(xm,x1,y1)

    rv      = @views @. sqrt((xm[:,1]-x1)^2+(xm[:,2]-y1)^2)
    return  @. cos(1.5*pi*rv)^3 * (rv<=1/3)

end

#-------------------------------------------------------------------------------
#VELOCITY FIELDS:

a_Radial(t,x,K,omega,lambdax,lambday)   = @views hcat(K.*sin.(lambdax.*x[:,1]).*cos(omega.*t),
                                                    K.*sin.(lambday.*x[:,2]).*cos(omega.*t) )
a_Radial(t,x,Kx,Ky,omega,lambdax,lambday)= @views hcat(Kx.*sin.(lambdax.*x[:,1]).*cos(omega.*t),
                                                    Ky.*sin.(lambday.*x[:,2]).*cos(omega.*t) )
a_Radial_sin(t,x,Kx,Ky,omega,lambdax,lambday)= @views hcat(Kx.*sin.(lambdax.*x[:,1]).*sin(omega.*t),
                                                    Ky.*sin.(lambday.*x[:,2]).*sin(omega.*t) )
a_RigidBodyRotation(t,xm,omega)         = @views hcat(-omega*xm[:,2],omega*xm[:,1])

diva_Radial(t,x,Kx,Ky,omega,lambdax,lambday)    = @views @. Kx*lambdax*cos(lambdax*x[:,1])*cos(omega*t) +
                                                    Ky*lambday*cos(lambday*x[:,2])*cos(omega*t)

function a_TestIV(t,xm,epsilon,u0,omega,lambda)

    # Cxv         = @views @. -(epsilon*lambda*sin(omega*t)+omega/lambda*cos(omega*t))*sin(lambda*xm[:,1])
    # return hcat( Cxv./utheor_CD_testIV(t,xm,epsilon,u0,omega,lambda), zeros(size(xm,1)) )

    auxv        = -(epsilon*lambda*sin(omega*t)+omega/lambda*cos(omega*t))./utheor_CD_testIV(t,xm,epsilon,u0,omega,lambda)
    return      hcat( auxv.*sin.(lambda*xm[:,1]), auxv.*sin.(lambda*xm[:,2]) )

end

function a_TestIVb(t,xm,epsilon,u0,omega,lambda)

    # Cxv         = @views @. -(epsilon*lambda*sin(omega*t)+omega/lambda*cos(omega*t))*sin(lambda*xm[:,1])
    # return hcat( Cxv./utheor_CD_testIV(t,xm,epsilon,u0,omega,lambda), zeros(size(xm,1)) )

    auxv        = (epsilon*lambda*sin(omega*t)+omega/lambda*cos(omega*t))./utheor_CD_testIVb(t,xm,epsilon,u0,omega,lambda)
    return      hcat( auxv.*sin.(lambda*xm[:,1]), auxv.*sin.(lambda*xm[:,2]) )

end

function diva_TestIV(t,xm,epsilon,u0,omega,lambda)

    Ev          = -(epsilon*lambda*sin(omega*t)+omega/lambda*cos(omega*t))
    uv          = utheor_CD_testIV(t,xm,epsilon,u0,omega,lambda)
    axxv        = @views @. (Ev*lambda*cos(lambda*xm[:,1])*uv - Ev*sin(lambda*xm[:,1])*(-lambda)*sin(omega*t)*sin(lambda*xm[:,1]))/uv^2
    ayyv        = @views @. (Ev*lambda*cos(lambda*xm[:,2])*uv - Ev*sin(lambda*xm[:,2])*(-lambda)*sin(omega*t)*sin(lambda*xm[:,2]))/uv^2
    return      axxv .+ ayyv

end

function diva_TestIVb(t,xm,epsilon,u0,omega,lambda)

    Ev          = epsilon*lambda*sin(omega*t)+omega/lambda*cos(omega*t)
    uv          = utheor_CD_testIVb(t,xm,epsilon,u0,omega,lambda)
    axxv        = @views @. (Ev*lambda*cos(lambda*xm[:,1])*uv - Ev*sin(lambda*xm[:,1])*(+lambda)*sin(omega*t)*sin(lambda*xm[:,1]))/uv^2
    ayyv        = @views @. (Ev*lambda*cos(lambda*xm[:,2])*uv - Ev*sin(lambda*xm[:,2])*(+lambda)*sin(omega*t)*sin(lambda*xm[:,2]))/uv^2
    return      axxv .+ ayyv

end

#-------------------------------------------------------------------------------
#INCOMING FLOW DATA:

uOut_null(t,xm)     = zeros(size(xm,1))

#-------------------------------------------------------------------------------
#BOUNDARY CONDITIONS:

function uq_testIII(boundary,t,xm,epsilon)
    Nx          = size(xm,1)
    if boundary==1
        return @views @. epsilon*(0.5+0.5*cos(pi*xm[:,1]))
    elseif boundary==2
        return zeros(Nx)
    elseif boundary==3
        return @views @. -epsilon*(0.5+0.5*cos(pi*xm[:,1]))
    elseif boundary==4
        return ones(Nx)
    end
end

function uq_testIIIb(boundary,t,xm,epsilon,lambda)
    Nx          = size(xm,1)
    if boundary==1
        return @views @. 2*real(1+exp(lambda*(xm[:,1]+1)));
    elseif boundary==2
        return @views @. -epsilon*real(lambda*exp(lambda*(1+1)))*real(1+exp(lambda*(xm[:,2]+1)));
    elseif boundary==3
        return @views @. -epsilon*real(1+exp(lambda*(xm[:,1]+1)))*real(lambda*exp(lambda*(1+1)));
    elseif boundary==4
        return @views @. 2*real(1+exp(lambda*(xm[:,2]+1)))
    end
end

function uq_testIIIc(boundary,t,xm,epsilon)
    Nx          = size(xm,1)
    if boundary==1
        # return @views @. epsilon*(0.5+0.5*cos(pi*xm[:,1]))
        return @views @. epsilon*exp(-xm[:,1]-1)
    elseif boundary==2
        # return zeros(Nx)
        return @views @. epsilon*exp(-2)*(2+xm[:,2])
    elseif boundary==3
        # return @views @. -epsilon*(0.5+0.5*cos(pi*xm[:,1]))
        return @views @. -epsilon*exp(-xm[:,1]-1)
    elseif boundary==4
        # return ones(Nx)
        return @views @. 2+xm[:,2]
    end
end

function uq_testIV(boundary,t,xm,epsilon,u0,omega,lambda)
    Nx          = size(xm,1)
    if boundary==1
        # return zeros(Nx)
        return @views @. u0 + sin(omega*t)*cos(lambda*xm[:,1]) + sin(omega*t)*cos(lambda*xm[:,2])
    elseif boundary==2
        return @views @. -epsilon*sin(omega*t)*lambda*(-sin(lambda*xm[:,1]))
    elseif boundary==3
        # return zeros(Nx)
        return @views @. -epsilon*sin(omega*t)*lambda*(-sin(lambda*xm[:,2]))
    elseif boundary==4
        # return ones(Nx)
        return @views @. u0 + sin(omega*t)*cos(lambda*xm[:,1]) + sin(omega*t)*cos(lambda*xm[:,2])
    end
end

function uq_testIVb(boundary,t,xm,epsilon,u0,omega,lambda)
    Nx          = size(xm,1)
    if boundary==1
        # return zeros(Nx)
        return @views @. u0 - sin(omega*t)*cos(lambda*xm[:,1]) - sin(omega*t)*cos(lambda*xm[:,2])
    elseif boundary==2
        return @views @. -epsilon*sin(omega*t)*lambda*sin(lambda*xm[:,1])
    elseif boundary==3
        # return zeros(Nx)
        return @views @. -epsilon*sin(omega*t)*lambda*sin(lambda*xm[:,2])
    elseif boundary==4
        # return ones(Nx)
        return @views @. u0 - sin(omega*t)*cos(lambda*xm[:,1]) - sin(omega*t)*cos(lambda*xm[:,2])
    end
end

#-------------------------------------------------------------------------------
#THEORETICAL SOLUTIONS:

function utheor_CD_testIV(t,xm,epsilon,u0,omega,lambda)

    # return @views @. u0 + sin(omega*t)*cos(lambda*xm[:,1])
    return @views @. u0 + sin(omega*t)*cos(lambda*xm[:,1]) + sin(omega*t)*cos(lambda*xm[:,2])

end

function utheor_CD_testIVb(t,xm,epsilon,u0,omega,lambda)

    # return @views @. u0 + sin(omega*t)*cos(lambda*xm[:,1])
    return @views @. u0 - sin(omega*t)*cos(lambda*xm[:,1]) - sin(omega*t)*cos(lambda*xm[:,2])

end

function utheor_Radial(u0fun,t,xm,K,omega,lambdax,lambday)

    x0m                     = zeros(size(xm,1),2)
    @views @. x0m[:,1]      = sign(xm[:,1])/lambdax*acos(tanh(atanh(cos(lambdax*xm[:,1]))+lambdax*K/omega*sin(omega*t)))
    @views @. x0m[:,2]      = sign(xm[:,2])/lambday*acos(tanh(atanh(cos(lambday*xm[:,2]))+lambday*K/omega*sin(omega*t)))
    Jv                      = zeros(size(xm,1))
    @views @. Jv            = ( (abs(sin(lambdax*x0m[:,1]))+1e-10)/(abs(sin(lambdax*xm[:,1]))+1e-10) )*
                              ( (abs(sin(lambday*x0m[:,2]))+1e-10)/(abs(sin(lambday*xm[:,2]))+1e-10) )
    u0v                     = u0fun(x0m)
    uv                      = u0v .* Jv
    return                  uv

end

function utheor_Radial(u0fun,t,xm,Kx,Ky,omega,lambdax,lambday)

    x0m                     = zeros(size(xm,1),2)
    @views @. x0m[:,1]      = sign(xm[:,1])/lambdax*acos(tanh(atanh(cos(lambdax*xm[:,1]))+lambdax*Kx/omega*sin(omega*t)))
    @views @. x0m[:,2]      = sign(xm[:,2])/lambday*acos(tanh(atanh(cos(lambday*xm[:,2]))+lambday*Ky/omega*sin(omega*t)))
    Jv                      = zeros(size(xm,1))
    @views @. Jv            = ( (abs(sin(lambdax*x0m[:,1]))+1e-10)/(abs(sin(lambdax*xm[:,1]))+1e-10) )*
                              ( (abs(sin(lambday*x0m[:,2]))+1e-10)/(abs(sin(lambday*xm[:,2]))+1e-10) )
    u0v                     = u0fun(x0m)
    uv                      = u0v .* Jv
    return                  uv

end

function utheor_Radial_sin(u0fun,t,xm,Kx,Ky,omega,lambdax,lambday)

    x0m                     = zeros(size(xm,1),2)
    @views @. x0m[:,1]      = sign(xm[:,1])/lambdax*acos(tanh(atanh(cos(lambdax*xm[:,1]))+lambdax*Kx/omega*(1-cos(omega*t))))
    @views @. x0m[:,2]      = sign(xm[:,2])/lambday*acos(tanh(atanh(cos(lambday*xm[:,2]))+lambday*Ky/omega*(1-cos(omega*t))))
    Jv                      = zeros(size(xm,1))
    @views @. Jv            = ( (abs(sin(lambdax*x0m[:,1]))+1e-10)/(abs(sin(lambdax*xm[:,1]))+1e-10) )*
                              ( (abs(sin(lambday*x0m[:,2]))+1e-10)/(abs(sin(lambday*xm[:,2]))+1e-10) )
    u0v                     = u0fun(x0m)
    uv                      = u0v .* Jv
    return                  uv

end



function utheor_RigidBodyRotation(u0fun,t,xm,omega)

    x0m                     = zeros(size(xm,1),2)
    @views @. x0m[:,1]      = @views @. xm[:,1]*cos(omega*t)+xm[:,2]*sin(omega*t)
    @views @. x0m[:,2]      = @views @. xm[:,2]*cos(omega*t)-xm[:,1]*sin(omega*t)
    u0v                     = u0fun(x0m)
    return                  u0v

end

function utheor_Burgers_RampHeaviside(t,xm,xramp,xHeav)
    x   = view(xm,:,1)
    u   = @. (x-xramp)/(1*t)*(x>xramp)*(x<=xramp+1*t) +
                    1*(x>xramp+1*t)*(x<xHeav+0.5*t) + 0.5*(x==xHeav+0.5*t);
    return u
end

function utheor_Burgers_DoubleRamp(t,xm,x1,x2,x3,x4,y1,y2,y3)

    x           = view(xm,:,1)
    u           = zeros(size(x))
    tshock      = (x4-x3)/(y2-y3)
    xshock      = x3+y2*tshock
    for ii=1:length(x)
        if x[ii]<=x1+y1*t
            u[ii]   = y1
        elseif x[ii]<=x2+y2*t
            u[ii]   = y1 + (y2-y1)/(x2+y2*t-x1-y1*t)*(x[ii]-(x1+y1*t))
        else
            if t<=tshock
                #No shock yet
                if x[ii]<=x3+y2*t
                    u[ii]   = y2
                elseif x[ii]<=x4+y3*t
                    u[ii]   = y2 + (y3-y2)/(x4+y3*t-x3-y2*t)*(x[ii]-(x3+y2*t))
                else
                    u[ii]   = y3
                end
            else
                #Shock:
                if x[ii]<=xshock+0.5*(y2+y3)*(t-tshock)
                    u[ii]   = y2
                else
                    u[ii]   = y3
                end
            end
        end
    end

    return u

end

function utheor_Burgers_Tanh(t,xm,A,delta,B)

    #Compute departure points:
    F(xi)       = A*tanh(xi/delta)+B
    Fdot(xi)    = A/delta/cosh(xi/delta)^2
    xiv         = zeros(size(xm,1))
    for ii=1:size(xm,1)
        x       = xm[ii,1]
        xi_nm1  = x
        iter    = 0
        while true
            xi_n    = xi_nm1 - (F(xi_nm1)*t+xi_nm1-x)/(Fdot(xi_nm1)*t+1)
            iter    += 1
            if norm(xi_n-xi_nm1)<1e-10
                xiv[ii]     = xi_n
                break
            elseif iter==100
                @warn "Not converged when computing feet of the characteristics"
                xiv[ii]     = xi_n
                break
            else
                xi_nm1      = xi_n
            end
        end
    end

    #Return u:
    return F.(xiv)

end

function utheor_CD_TestIII(xm)

    #Get matrix with coefficients:
    av      = readdlm("$(@__DIR__)/TestDirNeum_sol.txt")
    Nx      = 81
    Ny      = 81
    am      = reshape(av,Nx,Ny)

    #Loop points in xm:
    nPoints = size(xm,1)
    uv      = zeros(nPoints)
    Pv      = zeros(Nx)
    Qv      = zeros(Ny)
    for ii=1:nPoints

        #Compute Chebyshev polynomials:
        Pv[1]       = 1
        Pv[2]       = xm[ii,1]
        [ Pv[ix]    = 2*xm[ii,1]*Pv[ix-1] - Pv[ix-2] for ix=3:Nx ]
        Qv[1]       = 1
        Qv[2]       = xm[ii,2]
        [ Qv[ix]    = 2*xm[ii,2]*Qv[ix-1] - Qv[ix-2] for ix=3:Nx ]

        #Compute u:
        uv[ii]      = dot(Pv,am*Qv)

    end

    return uv

end

function utheor_CD_TestIIIb(xm)

    #Get matrix with coefficients:
    av      = readdlm("$(@__DIR__)/TestIIIb.txt")
    Nx      = Int(sqrt(length(av)))
    Ny      = Nx
    am      = reshape(av,Nx,Ny)

    #Loop points in xm:
    nPoints = size(xm,1)
    uv      = zeros(nPoints)
    Pv      = zeros(Nx)
    Qv      = zeros(Ny)
    for ii=1:nPoints

        #Compute Chebyshev polynomials:
        Pv[1]       = 1
        Pv[2]       = xm[ii,1]
        [ Pv[ix]    = 2*xm[ii,1]*Pv[ix-1] - Pv[ix-2] for ix=3:Nx ]
        Qv[1]       = 1
        Qv[2]       = xm[ii,2]
        [ Qv[ix]    = 2*xm[ii,2]*Qv[ix-1] - Qv[ix-2] for ix=3:Nx ]

        #Compute u:
        uv[ii]      = dot(Pv,am*Qv)

    end

    return uv

end

function utheor_CD_TestIIIc(xm)

    #Get matrix with coefficients:
    av      = readdlm("$(@__DIR__)/TestIIIc.txt")
    Nx      = Int(sqrt(length(av)))
    Ny      = Nx
    am      = reshape(av,Nx,Ny)

    #Loop points in xm:
    nPoints = size(xm,1)
    uv      = zeros(nPoints)
    Pv      = zeros(Nx)
    Qv      = zeros(Ny)
    for ii=1:nPoints

        #Compute Chebyshev polynomials:
        Pv[1]       = 1
        Pv[2]       = xm[ii,1]
        [ Pv[ix]    = 2*xm[ii,1]*Pv[ix-1] - Pv[ix-2] for ix=3:Nx ]
        Qv[1]       = 1
        Qv[2]       = xm[ii,2]
        [ Qv[ix]    = 2*xm[ii,2]*Qv[ix-1] - Qv[ix-2] for ix=3:Nx ]

        #Compute u:
        uv[ii]      = dot(Pv,am*Qv)

    end

    return uv

end

#-------------------------------------------------------------------------------
#NONLINEAR FLUX:

function isDir_Burgers2D_Rectangle(t::Float64,xm::Array{Float64,2},boundary::Int)
    Nx      = size(xm,1)
    if boundary==1 || boundary==4
        return ones(Bool,Nx)
    else
        return zeros(Bool,Nx)
    end
end

function isDir_RigidBody_Rectangle(t::Float64,xm::Array{Float64,2},boundary::Int)
    Nx      = size(xm,1)
    if boundary==1
        return @views @. xm[:,1]>0
    elseif boundary==2
        return @views @. xm[:,2]>0
    elseif boundary==3
        return @views @. xm[:,1]<0
    elseif boundary==4
        return @views @. xm[:,2]<0
    end
end
