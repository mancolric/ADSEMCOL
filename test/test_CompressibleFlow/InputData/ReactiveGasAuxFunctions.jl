function ChapmanJouget(gamma::Float64, q::Float64)

    #We assume the state at the end p2=1, rho2=1, u2=sqrt(gamma), YF2=0.
    p2      = 1.0
    rho2    = 1.0
    u2      = sqrt(gamma*p2/rho2)
    
    #Then, we compute the state "0" by imposing conservation equations:
    #   rho0*u0                                                 = M
    #   p0+rho0*u0^2                                            = P
    #   (gamma/(gamma-1)*p0 + 0.5*rho0*u0^2 + rho0*YF0*Q0)*u0   = HU
    #Note that rho0*Q0*YF0 = q*gamma/(gamma-1)*p0, so the last equation is
    #   ((1+q)*gamma/(gamma-1)*p0 + 0.5*rho0*u0^2)*u0           = HU
    #If we write the equations in terms of u0:
    #A u0^2 + B*u0 + HU = 0
    M   = rho2*u2
    P   = p2+rho2*u2^2
    HU  = (gamma/(gamma-1.0)*p2 + 0.5*rho2*u2^2)*u2
    #
    A   = ((1.0+q)*gamma/(gamma-1.0)-0.5)*M
    B   = -(1.0+q)*gamma/(gamma-1.0)*P
    C   = HU
    #
    u0      = (-B+sqrt(B^2-4*A*C))/(2*A)
    rho0    = M/u0
    p0      = P-M*u0
    M0      = u0/sqrt(gamma*p0/rho0)
    
    #Rescale solution with entry conditions:
    rho2    *= 1.0/rho0
    p2      *= 1.0/p0
    u2      = sqrt(gamma*p2/rho2)
    #
    rho0    = 1.0
    p0      = 1.0
    u0      = M0*sqrt(gamma)
    
    #=
    #Check results:
    display(rho0*u0 - rho2*u2)
    display(p0+rho0*u0^2 - (p2+rho2*u2^2))
    display(((1+q)*gamma/(gamma-1)*p0 + 0.5*rho0*u0^2)*u0 - 
            (gamma/(gamma-1)*p2 + 0.5*rho2*u2^2)*u2)
    =#
    
    #State at region 1:
    rho1            = rho0 / ( (gamma-1.0)/(gamma+1.0) + 2.0/(gamma+1.0)/M0^2 )
    u1              = rho0*u0/rho1
    p1              = p0+rho0*u0^2-rho1*u1^2
    
    return M0, rho0, u0, p0, rho1, u1, p1, rho2, u2, p2
    
end

function FP_ChapmanJouget(gamma::Float64, q::Float64, beta0::Float64; 
    Deltax::Float64=1e-3, xf::Float64=20.0, YF0::Float64=0.01, PlotRes::Bool=false)

    #Compute detonation states:
    M0, rho0, u0, p0, rho1, u1, p1, rho2, u2, p2 = ChapmanJouget(gamma,q)
    
    #Compute constants:
    RT0             = p0/rho0
    RT1             = p1/rho1
    RTa             = beta0*RT0                     #Activation temperature
    beta1           = RTa/RT1                       #beta coeff with respect to CJ-1 state
    Q               = gamma/(gamma-1)*RT0/YF0*q     #Fuel formation mass enthalpy 
    BB              = u1*exp(beta1)/(1.0*beta1)     #Characteristic length is 1
                        
    #Function to compute rho, u, p, YF from vector
    #   F=[rhoY_F u, rho u, p + rho u^2, ( gamma/(gamma-1)*p + rhoY_F Q + 1/2 rho u^2 ) u ]
    function ExtractVars(Fv::Vector{Float64})
    
        #YF = rho YF u / (rho u):
        YF          = Fv[1]/Fv[2]
        
        #F2 = rho u
        #F3 = p + rho u^2 
        #F4 = ( gamma/(gamma-1)*p + 1/2 rho u^2 + rho Q YF ) u
        #   = gamma/(gamma-1)*p*u + 1/2 F2 u^2 + F1 Q
        #Combining these equations:
        #(gamma/(gamma-1)-1/2) F2 U^2 - gamma/(gamma-1) F3 U + F4-F1*Q = 0
        A       = (gamma/(gamma-1)-0.5)*Fv[2]
        B       = -gamma/(gamma-1)*Fv[3]
        C       = Fv[4]-Fv[1]*Q
        u       = (-B-sqrt(B^2-4*A*C))/(2*A)    #Positive root provides state 0, negative provides 1
        
        #Therefore:
        rho         = Fv[2]/u
        p           = Fv[3]-Fv[2]*u
        
        return rho, u, p, YF
        
    end
    
    #Function to march flux variables in space:
    function dF_dx(Fv::Vector{Float64})
    
        rho, u, p, YF   = ExtractVars(Fv)
        RT              = p/rho
        mFdot           = -BB*rho*YF*exp(-RTa/RT)
        return [ mFdot; 0; 0; 0 ]
        
    end
    
    #Initialize variables:
    x_n             = 0.0
    F_n             = [ rho0*YF0*u0 ;
                        rho0*u0;
                        p0 + rho0*u0^2;
                        (gamma/(gamma-1.0)*p0 + rho0*YF0*Q + 0.5*rho0*u0^2)*u0 ]
                        
#     display((rho0, u0, p0, YF0))
#     display((rho1, u1, p1, YF0))
#     display(ExtractVars(F_n))
#     return
    
    #Apply RK method:
    RK              = RK_Coefficients("Ascher3")
    NF              = length(F_n)
    fm              = zeros(NF, RK.stages)
    F_np1           = zeros(NF)
    xv              = [ -xf*0.1, x_n, x_n ]
    rhov            = [ rho0, rho0, rho1 ]
    uv              = [ u0, u0, u1 ]
    pv              = [ p0, p0, p1 ]
    YFv             = [ YF0, YF0, YF0 ]
    while x_n < xf
    
        #Compute next stage:
        fm[:,1]         = dF_dx(F_n)
        for ii=2:RK.stages
            S           = zeros(NF)
            for jj=1:ii-1
                S       += RK.AE[ii,jj]*fm[:,jj]
            end
            F_np1       .= F_n + Deltax*S
            fm[:,ii]    = dF_dx(F_np1)
        end
        
        #Update solution:
        x_n         += Deltax
        F_n         .= F_np1
        
        #Save solution:
        xv                      = vcat(xv, x_n)
        rho_n, u_n, p_n, YF_n   = ExtractVars(F_n)
        rhov                    = vcat(rhov, rho_n)
        uv                      = vcat(uv, u_n)
        pv                      = vcat(pv, p_n)
        YFv                     = vcat(YFv, YF_n)
        
    end
    
    #Compute Mach number:
    Mv      = @. uv/sqrt(gamma*pv/rhov)
    
    #Plot results:
    if PlotRes
    
        Nx      = length(xv)
        display("theoretical")
        display((rho2, u2, p2, 0.0, 1.0))
        display("numerical")
        display((rhov[Nx], uv[Nx], pv[Nx], YFv[Nx], Mv[Nx]))
        
        figure()
        plot(xv, rhov, "b")
        plot(xv, uv, "g")
        plot(xv, pv, "c")
        plot(xv, pv./rhov, "r")
        plot(xv, YFv./YF0, color="orange")
        plot(xv, Mv, "k")
        legend([latexstring("\\rho"), 
                latexstring("u"), 
                latexstring("p"),
                latexstring("T"), 
                latexstring("Y_{F}/Y_{F0}"), 
                latexstring("M") ], 
                loc="best")
        xlabel("x/L") #Here, L=u1 exp(beta1)/(beta1*BB)
        grid("on")
        
    end
    
#     figure()
#     plot(xv, rhov.*uv.*YFv, "b")
#     plot(xv, rhov.*uv, "g")
#     plot(xv, pv+rhov.*uv.^2, "r")
#     plot(xv, gamma/(gamma-1)*pv.*uv + 0.5*rhov.*uv.^3 + rhov.*YFv.*Q.*uv, "k")
    
    return xv, rhov, uv, pv, YFv
    
end

function FXP_ChapmanJouget(; gamma::Float64=1.4, q::Float64=5.0, 
    RThatI::Float64=10.0, RThatB::Float64=1.0, #RThat=RT_a/RT_1, #1: means post-shock condition
    BIhat::Float64=2.5e-3, BBhat::Float64=4.4e5, 
    Deltax::Float64=1e-3, xf::Float64=20.0, YF0::Float64=0.01, YX0::Float64=0.0, 
    PlotRes::Bool=false)

    #Variables will be make dimensionless via p1, rho1 and B_R:
    
    #Compute detonation states:
    M0, rho0, u0, p0, rho1, u1, p1, rho2, u2, p2 = ChapmanJouget(gamma,q)
    
    #Characteristic variables:
    rhoc            = rho1
    uc              = u1
    pc              = rhoc*uc^2
    #tc=1 via BR
    
    #Scale pre and post shock conditions:
#     rho0            /= rhoc
#     u0              /= uc
#     p0              /= pc
#     rho1            /= rhoc
#     u1              /= uc
#     p1              /= pc
#     rho2            /= rhoc
#     u2              /= uc
#     p2              /= pc
    
    #Compute RT:
    RT0             = p0/rho0
    RT1             = p1/rho1
    
    #Compute constants:
    Q               = gamma/(gamma-1)*RT0/YF0*q     #F and X formation mass enthalpy 
    RTI             = RThatI*RT1                    #Activation temperature
    RTB             = RThatB*RT1                    #Activation temperature
    BR              = 1.0                           #t_c = B_R^{-1} = 1.0
    BI              = BIhat * BR
    BB              = BBhat * BR/rho1
                        
    #Function to compute rho, u, p, YF, YX from vector 
    #   F   =   [   rhoY_F u, 
    #               rhoY_X u, 
    #               rho u, 
    #               p + rho u^2, 
    #               ( gamma/(gamma-1)*p - rho (Y_F+Y_X) Q + 1/2 rho u^2 ) u ]
    function ExtractVars(Fv::Vector{Float64})
    
        #YF = rho YF u / (rho u):
        YF          = Fv[1]/Fv[3]
        YX          = Fv[2]/Fv[3]
        YP          = 1-YF-YX
        
        #F3 = rho u
        #F4 = p + rho u^2 
        #F5 = ( gamma/(gamma-1)*p + 1/2 rho u^2 + rho Q (YF+YX) ) u
        #   = gamma/(gamma-1)*p*u + 1/2 F3 u^2 + (F1+F2) Q
        #Combining these equations:
        #(gamma/(gamma-1)-1/2) F3 U^2 - gamma/(gamma-1) F4 U + F5 - (F1+F2)*Q = 0
        A       = (gamma/(gamma-1)-0.5)*Fv[3]
        B       = -gamma/(gamma-1)*Fv[4]
        C       = Fv[5]-(Fv[1]+Fv[2])*Q
        u       = (-B-sqrt(B^2-4*A*C))/(2*A)    #Positive root provides state 0, negative provides 1
        
        #Therefore:
        rho         = Fv[3]/u
        p           = Fv[4]-Fv[3]*u
        
        return rho, u, p, YF, YX
        
    end
    
    #Function to march flux variables in space:
    function dF_dx(Fv::Vector{Float64})
    
        rho, u, p, YF, YX   = ExtractVars(Fv)
        RT                  = p/rho
        kI                  = BI*exp(-RTI/RT)
        kB                  = BB*exp(-RTB/RT)
        kR                  = BR
        mFdot               = - kI*rho*YF - kB*rho^2*YF*YX
        mXdot               =   kI*rho*YF + kB*rho^2*YF*YX - kR*rho*YX
        return [ mFdot; mXdot; 0; 0; 0 ]
        
    end
    
    #Initialize variables:
    x_n             = 0.0
    F_n             = [ rho0*YF0*u0 ;
                        rho0*YX0*u0 ;
                        rho0*u0;
                        p0 + rho0*u0^2;
                        (gamma/(gamma-1.0)*p0 + rho0*(YF0+YX0)*Q + 0.5*rho0*u0^2)*u0 ]
                        
#     display((rho0, u0, p0, YF0, YX0))
#     display((rho1, u1, p1, YF0, YX0))
#     display(ExtractVars(F_n))
#     return
    
    #Apply RK method:
    RK              = RK_Coefficients("Ascher3")
    NF              = length(F_n)
    fm              = zeros(NF, RK.stages)
    F_np1           = zeros(NF)
    xv              = [ -xf*0.1, x_n, x_n ]
    rhov            = [ rho0, rho0, rho1 ]
    uv              = [ u0, u0, u1 ]
    pv              = [ p0, p0, p1 ]
    YFv             = [ YF0, YF0, YF0 ]
    YXv             = [ YX0, YX0, YX0 ]
        
    while x_n < xf
    
        println("x_n=", x_n)
        
        #Compute next stage:
        fm[:,1]         = dF_dx(F_n)
        for ii=2:RK.stages
            S           = zeros(NF)
            for jj=1:ii-1
                S       += RK.AE[ii,jj]*fm[:,jj]
            end
            F_np1       .= F_n + Deltax*S
            fm[:,ii]    = dF_dx(F_np1)
        end
        
        #Update solution:
        x_n             += Deltax
        F_n             .= F_np1
        
        #Save solution:
        xv                              = vcat(xv, x_n)
        rho_n, u_n, p_n, YF_n, YX_n     = ExtractVars(F_n)
        rhov                            = vcat(rhov, rho_n)
        uv                              = vcat(uv, u_n)
        pv                              = vcat(pv, p_n)
        YFv                             = vcat(YFv, YF_n)
        YXv                             = vcat(YXv, YX_n)
        
    end
    
    #Compute Mach number:
    Mv      = @. uv/sqrt(gamma*pv/rhov)
    
    #Plot results:
    if false
    
        Nx      = length(xv)
        display("theoretical")
        display((rho2, u2, p2, 0.0, 1.0))
        display("numerical")
        display((rhov[Nx], uv[Nx], pv[Nx], YFv[Nx], Mv[Nx]))
        
        figure()
        plot(xv, rhov, "b")
        plot(xv, uv, "g")
        plot(xv, pv, "c")
        plot(xv, pv./rhov, "r")
        plot(xv, YFv, color="orange")
        plot(xv, YXv, color="yellow")
        plot(xv, Mv, "k")
        legend([latexstring("\\rho"), 
                latexstring("u"), 
                latexstring("p"),
                latexstring("T"), 
                latexstring("Y_{F}"), 
                latexstring("Y_{X}"),
                latexstring("M") ], 
                loc="best")
        xlabel("x/L") #Here, L=u1 exp(beta1)/(beta1*BB)
        grid("on")
        
    end
    #Plot each variable scaled with post-shock conditions (1):
    if PlotRes
    
        Nx      = length(xv)
        display("theoretical")
        display((rho2, u2, p2, 0.0, 1.0))
        display("numerical")
        display((rhov[Nx], uv[Nx], pv[Nx], YFv[Nx], Mv[Nx]))
        
        figure()
        plot(xv, rhov/rho1, "b")
        plot(xv, uv/u1, "g")
        plot(xv, pv/p1, "c")
        plot(xv, (pv./rhov)/RT1, "r")
        plot(xv, Mv, "k")
        legend([latexstring("\\rho"), 
                latexstring("u"), 
                latexstring("p"),
                latexstring("T"), 
                latexstring("M") ], 
                loc="best")
        xlabel("x/L") #Here, L=u1 exp(beta1)/(beta1*BB)
        grid("on")
        
        figure()
        plot(xv, YFv, color="b")
        plot(xv, YXv, color="orange")
        plot(xv, 1.0.-YFv.-YXv, color="g")
        legend([latexstring("Y_{F}"), 
                latexstring("Y_{X}"),
                latexstring("Y_{P}")], 
                loc="best")
        xlabel("x/L") #Here, L=u1 exp(beta1)/(beta1*BB)
        grid("on")
    
    end
    
    return xv, rhov, uv, pv, YFv, YXv
    
end
