using Interpolations

#-------------------------------------------------------------------------------
#INITIAL CONDITIONS:

function vx_SmoothVortex(xm::Array{Float64,2},gamma::Float64,beta::Float64,
    xvortex::Float64,yvortex::Float64,uinfty::Float64,vinfty::Float64)

    r2v     = @views @. (xm[:,1]-xvortex)^2+(xm[:,2]-yvortex)^2
    vx      = @views @. uinfty-(xm[:,2]-yvortex)*beta/2/pi*exp((1-r2v)/2)
    return  vx

end
function vy_SmoothVortex(xm::Array{Float64,2},gamma::Float64,beta::Float64,
    xvortex::Float64,yvortex::Float64,uinfty::Float64,vinfty::Float64)

    r2v     = @views @. (xm[:,1]-xvortex)^2+(xm[:,2]-yvortex)^2
    vy      = @views @. vinfty+(xm[:,1]-xvortex)*beta/2/pi*exp((1-r2v)/2)
    return  vy

end
function TStar_SmoothVortex(xm::Array{Float64,2},gamma::Float64,beta::Float64,
    xvortex::Float64,yvortex::Float64,uinfty::Float64,vinfty::Float64)

    r2v     = @views @. (xm[:,1]-xvortex)^2+(xm[:,2]-yvortex)^2
    Tv      = @views @. 1-(gamma-1)*beta^2/8/gamma/pi^2*exp(1-r2v)
    return  Tv

end
function e_SmoothVortex(xm::Array{Float64,2},gamma::Float64,beta::Float64,
    xvortex::Float64,yvortex::Float64,uinfty::Float64,vinfty::Float64)

    #e = cv*T = cv * Tstar * (p_infty/Rg/rho_infty) = 1/(gamma-1)*p_infty/rho_infty*Tstar
    return  TStar_SmoothVortex(xm,gamma,beta,xvortex,yvortex,uinfty,vinfty)./(gamma-1)

end
function rho_SmoothVortex(xm::Array{Float64,2},gamma::Float64,beta::Float64,
    xvortex::Float64,yvortex::Float64,uinfty::Float64,vinfty::Float64)

    return  TStar_SmoothVortex(xm,gamma,beta,xvortex,yvortex,uinfty,vinfty).^(1/(gamma-1))

end

function SodTube_Ms(rhoL::Float64,pL::Float64,rhoR::Float64,pR::Float64,gamma::Float64)

    aL      = sqrt(gamma*pL/rhoL)
    aR      = sqrt(gamma*pR/rhoR)
    fun     = (M) -> M-1/M-(gamma+1)/(gamma-1)*aL/aR*(1-(pR/pL*
                        (2*gamma/(gamma+1)*M^2-(gamma-1)/(gamma+1)))^((gamma-1)/2/gamma))
    dfun    = (M) -> 1+1/M^2+aL/aR*(pR/pL)^((gamma-1)/2/gamma)*(2*gamma/(gamma+1)*M^2-(gamma-1)/(gamma+1))^(-(gamma+1)/2/gamma)

    Mn      = 1.0
    iter    = 0
    while true

        Mnp1    = Mn - fun(Mn)/dfun(Mn)
        iter    += 1
        if abs(Mnp1-Mn) < 1e-6
            return Mnp1
        end
        if iter>=100
            @warn "Reached maximum of iterations"
            return Mnp1
        end
        Mn      = Mnp1

    end

end

function SodTube_rho(xm_t::Array{Float64,1},rhoL::Float64,pL::Float64,rhoR::Float64,pR::Float64,gamma::Float64,Ms::Float64)

    #Allocate output:
    rhov        = zeros(length(xm_t))

    #Get variables:
    aR          = sqrt(gamma*pR/rhoR)
    u1          = 2*aR/(gamma+1)*(Ms-1/Ms)
    p1          = pR*(2*gamma/(gamma+1)*Ms^2-(gamma-1)/(gamma+1))
    rho1        = rhoR/(2/(gamma+1)/Ms^2+(gamma-1)/(gamma+1))
    a1          = sqrt(gamma*p1/rho1)
    u2          = u1
    p2          = p1
    rho2        = rhoL*(p2/pL)^(1/gamma)
    a2          = sqrt(gamma*p2/rho2)
    aL          = sqrt(gamma*pL/rhoL)
    ushock      = Ms*aR

    #Left region:
    aux         = @. xm_t<-aL
    @.rhov[aux] = rhoL

    #Expansion fan:
    aux         = @. -aL<=xm_t<u2-a2
    a_aux       = @. (2*aL-(gamma-1)*xm_t[aux])/(gamma+1)
    @.rhov[aux] = rhoL*(a_aux/aL)^(2/(gamma-1))

    #Region 2:
    aux         = @. u2-a2<=xm_t<u2
    @.rhov[aux] = rho2

    #Region 1:
    aux         = @. u1<=xm_t<ushock
    @.rhov[aux] = rho1

    #Right region:
    aux         = @. ushock<=xm_t
    @.rhov[aux] = rhoR

    #Exit:
    return rhov

end

function SodTube_p(xm_t::Array{Float64,1},rhoL::Float64,pL::Float64,rhoR::Float64,pR::Float64,gamma::Float64,Ms::Float64)

    #Allocate output:
    pv        = zeros(length(xm_t))

    #Get variables:
    aR          = sqrt(gamma*pR/rhoR)
    u1          = 2*aR/(gamma+1)*(Ms-1/Ms)
    p1          = pR*(2*gamma/(gamma+1)*Ms^2-(gamma-1)/(gamma+1))
    rho1        = rhoR/(2/(gamma+1)/Ms^2+(gamma-1)/(gamma+1))
    a1          = sqrt(gamma*p1/rho1)
    u2          = u1
    p2          = p1
    rho2        = rhoL*(p2/pL)^(1/gamma)
    a2          = sqrt(gamma*p2/rho2)
    aL          = sqrt(gamma*pL/rhoL)
    ushock      = Ms*aR

    #Left region:
    aux         = @. xm_t<-aL
    @.pv[aux]   = pL

    #Expansion fan:
    aux         = @. -aL<=xm_t<u2-a2
    a_aux       = @. (2*aL-(gamma-1)*xm_t[aux])/(gamma+1)
    @.pv[aux]   = pL*(a_aux/aL)^(2*gamma/(gamma-1))

    #Region 2:
    aux         = @. u2-a2<=xm_t<u2
    @.pv[aux]   = p2

    #Region 1:
    aux         = @. u1<=xm_t<ushock
    @.pv[aux]   = p1

    #Right region:
    aux         = @. ushock<=xm_t
    @.pv[aux]   = pR

    #Exit:
    return pv

end

function SodTube_u(xm_t::Array{Float64,1},rhoL::Float64,pL::Float64,rhoR::Float64,pR::Float64,gamma::Float64,Ms::Float64)

    #Allocate output:
    uv        = zeros(length(xm_t))

    #Get variables:
    aR          = sqrt(gamma*pR/rhoR)
    u1          = 2*aR/(gamma+1)*(Ms-1/Ms)
    p1          = pR*(2*gamma/(gamma+1)*Ms^2-(gamma-1)/(gamma+1))
    rho1        = rhoR/(2/(gamma+1)/Ms^2+(gamma-1)/(gamma+1))
    a1          = sqrt(gamma*p1/rho1)
    u2          = u1
    p2          = p1
    rho2        = rhoL*(p2/pL)^(1/gamma)
    a2          = sqrt(gamma*p2/rho2)
    aL          = sqrt(gamma*pL/rhoL)
    ushock      = Ms*aR

    #Left region:
    aux         = @. xm_t<-aL
    @.uv[aux]   = 0

    #Expansion fan:
    aux         = @. -aL<=xm_t<u2-a2
    a_aux       = @. (2*aL-(gamma-1)*xm_t[aux])/(gamma+1)
    @.uv[aux]   = 2/(gamma-1)*(aL-a_aux)

    #Region 2:
    aux         = @. u2-a2<=xm_t<u2
    @.uv[aux]   = u2

    #Region 1:
    aux         = @. u1<=xm_t<ushock
    @.uv[aux]   = u1

    #Right region:
    aux         = @. ushock<=xm_t
    @.uv[aux]   = 0

    #Exit:
    return uv

end

function SteadySW_u(xv::Array{Float64,1},gamma::Float64,mu::Float64,
    rhoL::Float64,uL::Float64,pL::Float64,rhoR::Float64,uR::Float64,pR::Float64)

    #Constants:
    m           = rhoL*uL
    P           = pL + rhoL*uL^2
    H           = pL*uL + uL*(pL/(gamma-1.0) + 0.5*rhoL*uL^2)
    @assert abs(rhoR*uR-m)<1e-4 "m is not constant"
    @assert abs(pR + rhoR*uR^2 - P)<1e-4 "P is not constant"
    @assert abs(pR*uR + uR*(pR/(gamma-1.0) + 0.5*rhoR*uR^2) - H)<1e-4 "H is not constant"
    mubar       = 4.0*mu/3.0/m
    delta       = 0.5*(gamma-1.0)

    #Allocate output:
    omegav      = zeros(size(xv))

    #Loop and find omega with bisection method:
    omegaL      = m*uL/P
    omegaR      = m*uR/P
    omega0      = max(omegaL,omegaR)
    omega1      = min(omegaL,omegaR)
    B1          = omega0/(omega0-omega1)
    C1          = omega1/(omega0-omega1)
    x0          = mubar/(1+delta)*(C1-B1)*log((omega0-omega1)/2)
    MaxIter     = 100
    for ii=1:length(xv)
        X               = (1.0+delta)/mubar*(xv[ii]-x0)
        omegaA          = omega1
        omegaB          = omega0
        omega           = NaN
        resA            = B1*log(omega0-omegaA) - C1*log(omegaA-omega1) - X
        resB            = B1*log(omega0-omegaB) - C1*log(omegaB-omega1) - X
        nIter           = 0
        while true
            #New omega:
            omega       = 0.5*(omegaA+omegaB)
            res         = B1*log(omega0-omega) - C1*log(omega-omega1) - X
            nIter       += 1
            #Exit loop:
            if abs(omegaB-omegaA)<1e-5
                break
            elseif nIter==MaxIter
                @warn "Unable to achieve convergence"
                break
            end
            #Update A and B:
            if resA*res<0
                omegaB  = omega
                resB    = res
            else
                omegaA  = omega
                resA    = res
            end
        end
        omegav[ii]      = omega
    end

    return P*omegav/m

end

function SteadySW_rho(xv::Array{Float64,1},gamma::Float64,mu::Float64,
    rhoL::Float64,uL::Float64,pL::Float64,rhoR::Float64,uR::Float64,pR::Float64)

    uv  = SteadySW_u(xv,gamma,mu,rhoL,uL,pL,rhoR,uR,pR)
    return (rhoL*uL)./uv

end

function SteadySW_p(xv::Array{Float64,1},gamma::Float64,mu::Float64,
    rhoL::Float64,uL::Float64,pL::Float64,rhoR::Float64,uR::Float64,pR::Float64)

    #Constants:
    m           = rhoL*uL
    P           = pL + rhoL*uL^2
    H           = pL*uL + uL*(pL/(gamma-1.0) + 0.5*rhoL*uL^2)
    @assert abs(rhoR*uR-m)<1e-4 "m is not constant"
    @assert abs(pR + rhoR*uR^2 - P)<1e-4 "P is not constant"
    @assert abs(pR*uR + uR*(pR/(gamma-1.0) + 0.5*rhoR*uR^2) - H)<1e-4 "H is not constant"
    mubar       = 4.0*mu/3.0/m
    delta       = 0.5*(gamma-1.0)

    #Get omega and thetav:
    omegav      = SteadySW_u(xv,gamma,mu,rhoL,uL,pL,rhoR,uR,pR)*m/P
    thetav      = @. delta*(omegav^2-2*omegav+2*H*m/P^2)

    #Compute pressures: p = rho*R*T = theta*P/omega
    return @. thetav*P/omegav

    #Compute pressures (other way). Note that
    # p+rho*u^2-4/3*mu*du/dx=P,
    # p = P - m*u - 4/3*mu*du/dx = P - m*omega*P/m - 4/3*mu*P/m*domega/dx =
    #   = P - omega*P - mubar*P*domega/dx.
    #From B1*ln(omega0-omega)-C1*ln(omega-omega1)=(1+delta)/mubar*x, we obtain
    # (-B1/(omega0-omega)-C1/(omega-omega1))*domega/dx = (1+delta)/mubar
    #and then
    # p = P -omega*P - P*(1+delta)/(-B1/(omega0-omega)-C1/(omega-omega1))
    omegaL      = m*uL/P
    omegaR      = m*uR/P
    omega0      = max(omegaL,omegaR)
    omega1      = min(omegaL,omegaR)
    B1          = omega0/(omega0-omega1)
    C1          = omega1/(omega0-omega1)
    pv          = @. P - omegav*P + P*(1.0+delta)/(-B1/(omega0-omegav)-C1/(omegav-omega1))

end

#frho * rho1, where rho1 is the post-shock density, provides the theoretical solution:
function Sedov2D_frho(lambdav::Array{Float64,1})

    auxm     = [1.0     1.0     1.0     1.0
                0.9998	0.9996	0.9972	0.9984
                0.9802	0.9645	0.7651	0.8658
                0.9644	0.9374	0.6281	0.7829
                0.9476	0.9097	0.5161	0.7122
                0.9295	0.8812	0.4233	0.6513
                0.9096	0.8514	0.345	0.5982
                0.8725	0.7999	0.2427	0.5266
                0.8442	0.7638	0.1892	0.4884
                0.8094	0.7226	0.1415	0.4545
                0.7629	0.672	0.0974	0.4241
                0.7242	0.6327	0.0718	0.4074
                0.6894	0.599	0.0545	0.3969
                0.639	0.5521	0.0362	0.3867
                0.5745	0.4943	0.0208	0.3794
                0.518	0.4448	0.0123	0.376
                0.4748	0.4074	0.0079	0.3746
                0.4222	0.362	0.0044	0.3737
                0.3654	0.3133	0.0021	0.3732
                0.3	    0.2572	0.0008	0.373
                0.25	0.2143	0.0003	0.3729
                0.2	    0.1714	0.0001	0.3729
                0.15	0.1286	0.0	    0.3729
                0.1	    0.0857	0.0	    0.3729
                0.05    0.05    0.0     0.3729
                0.0     0.0     0.0     0.3729 ]

    gamma           = 1.4
    lambda_terp     = view(auxm, :, 1)
    f_terp          = view(auxm, :, 3)
    fv              = zeros(length(lambdav))
    for ii=1:length(lambdav)
        lambda      = lambdav[ii]
        if lambda==0.0
            fv[ii]  = 0.0
        elseif lambda<1.0
            i1      = findlast(lambda_terp.>=lambda)
            i2      = i1+1
            fv[ii]  = f_terp[i1]*(lambda-lambda_terp[i2])/(lambda_terp[i1]-lambda_terp[i2]) + 
                      f_terp[i2]*(lambda-lambda_terp[i1])/(lambda_terp[i2]-lambda_terp[i1])
        else
            fv[ii]  = (gamma-1.0)/(gamma+1.0)
        end
        # println("λ=$(lambdav[ii]); f=$(fv[ii])")
    end

    return fv

end

function Sedov2D_fv(lambdav::Array{Float64,1})

    auxm     = [1.0     1.0     1.0     1.0
                0.9998	0.9996	0.9972	0.9984
                0.9802	0.9645	0.7651	0.8658
                0.9644	0.9374	0.6281	0.7829
                0.9476	0.9097	0.5161	0.7122
                0.9295	0.8812	0.4233	0.6513
                0.9096	0.8514	0.345	0.5982
                0.8725	0.7999	0.2427	0.5266
                0.8442	0.7638	0.1892	0.4884
                0.8094	0.7226	0.1415	0.4545
                0.7629	0.672	0.0974	0.4241
                0.7242	0.6327	0.0718	0.4074
                0.6894	0.599	0.0545	0.3969
                0.639	0.5521	0.0362	0.3867
                0.5745	0.4943	0.0208	0.3794
                0.518	0.4448	0.0123	0.376
                0.4748	0.4074	0.0079	0.3746
                0.4222	0.362	0.0044	0.3737
                0.3654	0.3133	0.0021	0.3732
                0.3	    0.2572	0.0008	0.373
                0.25	0.2143	0.0003	0.3729
                0.2	    0.1714	0.0001	0.3729
                0.15	0.1286	0.0	    0.3729
                0.1	    0.0857	0.0	    0.3729
                0.05    0.05    0.0     0.3729
                0.0     0.0     0.0     0.3729 ]

    lambda_terp     = view(auxm, :, 1)
    f_terp          = view(auxm, :, 2)
    fv              = zeros(length(lambdav))
    for ii=1:length(lambdav)
        lambda      = lambdav[ii]
        if lambda==0.0
            fv[ii]  = 0.0
        elseif lambda<1.0
            i1      = findlast(lambda_terp.>=lambda)
            i2      = i1+1
            fv[ii]  = f_terp[i1]*(lambda-lambda_terp[i2])/(lambda_terp[i1]-lambda_terp[i2]) + 
                      f_terp[i2]*(lambda-lambda_terp[i1])/(lambda_terp[i2]-lambda_terp[i1])
        else
            fv[ii]  = 0.0
        end
        # println("λ=$(lambdav[ii]); f=$(fv[ii])")
    end

    return fv

end

function Sedov2D_fp(lambdav::Array{Float64,1})

    auxm     = [1.0     1.0     1.0     1.0
                0.9998	0.9996	0.9972	0.9984
                0.9802	0.9645	0.7651	0.8658
                0.9644	0.9374	0.6281	0.7829
                0.9476	0.9097	0.5161	0.7122
                0.9295	0.8812	0.4233	0.6513
                0.9096	0.8514	0.345	0.5982
                0.8725	0.7999	0.2427	0.5266
                0.8442	0.7638	0.1892	0.4884
                0.8094	0.7226	0.1415	0.4545
                0.7629	0.672	0.0974	0.4241
                0.7242	0.6327	0.0718	0.4074
                0.6894	0.599	0.0545	0.3969
                0.639	0.5521	0.0362	0.3867
                0.5745	0.4943	0.0208	0.3794
                0.518	0.4448	0.0123	0.376
                0.4748	0.4074	0.0079	0.3746
                0.4222	0.362	0.0044	0.3737
                0.3654	0.3133	0.0021	0.3732
                0.3	    0.2572	0.0008	0.373
                0.25	0.2143	0.0003	0.3729
                0.2	    0.1714	0.0001	0.3729
                0.15	0.1286	0.0	    0.3729
                0.1	    0.0857	0.0	    0.3729
                0.05    0.05    0.0     0.3729
                0.0     0.0     0.0     0.3729 ]

    lambda_terp     = view(auxm, :, 1)
    f_terp          = view(auxm, :, 4)
    fv              = zeros(length(lambdav))
    for ii=1:length(lambdav)
        lambda      = lambdav[ii]
        if lambda==0.0
            fv[ii]  = 0.0
        elseif lambda<1.0
            i1      = findlast(lambda_terp.>=lambda)
            i2      = i1+1
            fv[ii]  = f_terp[i1]*(lambda-lambda_terp[i2])/(lambda_terp[i1]-lambda_terp[i2]) + 
                      f_terp[i2]*(lambda-lambda_terp[i1])/(lambda_terp[i2]-lambda_terp[i1])
        else
            fv[ii]  = 0.0
        end
        # println("λ=$(lambdav[ii]); f=$(fv[ii])")
    end

    return fv

end

function ShuOsher_rho(xm::Array{Float64,2})

    auxm    = [ -5.001  3.854
                -4.982	3.854
                -2.925	3.854
                -2.807	3.838
                -2.712	3.781
                -2.653	3.748
                -2.583	3.707
                -2.524	3.977
                -2.476	4.059
                -2.394	4.043
                -2.335	4.018
                -2.276	3.993
                -2.229	3.969
                -2.123	3.903
                -1.969	3.805
                -1.851	3.74
                -1.71	3.658
                -1.568	3.617
                -1.533	3.895
                -1.498	4.051
                -1.403	4.067
                -1.297	4.018
                -1.156	3.936
                -0.979	3.83
                -0.825	3.699
                -0.613	3.617
                -0.542	3.74
                -0.519	3.871
                -0.483	3.993
                -0.377	4.067
                -0.212	4.01
                -0.071	3.895
                0.083	3.773
                0.2	3.691
                0.295	3.617
                0.413	3.633
                0.46	3.822
                0.495	3.977
                0.613	4.059
                0.696	3.977
                0.731	3.813
                0.755	3.674
                0.79	4.1
                0.814	4.313
                0.861	4.591
                0.92	4.288
                0.967	3.895
                0.991	3.519
                1.061	3.159
                1.097	3.486
                1.12	3.854
                1.203	4.231
                1.262	3.936
                1.321	3.445
                1.403	3.011
                1.439	3.38
                1.474	3.838
                1.498	4.173
                1.545	4.591
                1.592	4.329
                1.616	4.01
                1.651	3.707
                1.71	3.388
                1.769	3.895
                1.792	4.231
                1.863	4.64
                1.91	4.386
                1.934	4.043
                1.958	3.756
                1.981	3.486
                2.028	3.167
                2.087	3.478
                2.111	3.813
                2.146	4.051
                2.182	4.223
                2.264	3.92
                2.311	3.314
                2.358	3.036
                2.417	1.252
                2.441	0.908
                2.618	1.08
                2.736	1.154
                2.854	1.195
                2.972	1.154
                3.066	1.088
                3.16	1.007
                3.267	0.908
                3.361	0.818
                3.502	0.786
                3.632	0.827
                3.738	0.925
                3.821	1.015
                3.915	1.097
                3.986	1.154
                4.104	1.195
                4.233	1.146
                4.328	1.088
                4.399	1.007
                4.517	0.908
                4.611	0.835
                4.705	0.794
                4.811	0.81
                4.882	0.859
                4.982	0.908
                5.001   0.946 ]

    xv      = auxm[:,1]
    rhov    = auxm[:,2]

    #Interpolate:
    interp  = LinearInterpolation(xv, rhov)
    return interp(xm[:,1])

end

function ShuOsher_v(xm::Array{Float64,2})

    auxm    = [ -5.009	2.639
    -4.625	2.632
    -4.283	2.632
    -3.891	2.639
    -3.481	2.632
    -3.055	2.624
    -2.816	2.662
    -2.645	2.707
    -2.577	2.624
    -2.491	2.542
    -2.253	2.594
    -2.031	2.662
    -1.826	2.722
    -1.638	2.767
    -1.57	2.654
    -1.485	2.542
    -1.229	2.579
    -0.99	2.654
    -0.802	2.73
    -0.631	2.767
    -0.546	2.654
    -0.495	2.579
    -0.375	2.542
    -0.137	2.617
    0.085	2.692
    0.239	2.76
    0.392	2.73
    0.444	2.632
    0.529	2.564
    0.666	2.534
    0.904	2.609
    1.075	2.692
    1.246	2.76
    1.451	2.685
    1.553	2.579
    1.741	2.534
    1.928	2.602
    2.082	2.7
    2.184	2.767
    2.355	2.73
    2.44	0.018
    2.747	0.018
    3.02	0.018
    3.294	0.01
    3.584	0.01
    3.891	0.01
    4.232	0.01
    4.573	0.01
    5.009	0.01 ]

    xv      = auxm[:,1]
    vv      = auxm[:,2]

    #Interpolate:
    interp  = LinearInterpolation(xv, vv)
    return interp(xm[:,1])

end

function ShuOsher_p(xm::Array{Float64,2})

    auxm    = [ -5.001	10.298
    -4.428	10.298
    -3.908	10.298
    -3.44	10.298
    -3.007	10.326
    -2.799	10.128
    -2.574	9.702
    -2.539	11.121
    -2.314	10.752
    -2.106	10.326
    -1.828	9.759
    -1.586	9.277
    -1.516	11.092
    -1.17	10.582
    -0.962	9.986
    -0.615	9.305
    -0.546	10.099
    -0.494	11.035
    -0.39	11.092
    -0.182	10.752
    -0.009	10.27
    0.113	9.787
    0.321	9.333
    0.442	9.986
    0.511	10.638
    0.581	11.007
    0.685	11.064
    0.823	10.809
    0.945	10.44
    1.049	9.986
    1.153	9.617
    1.308	9.333
    1.412	9.929
    1.516	10.638
    1.603	11.007
    1.707	11.121
    1.846	10.894
    1.932	10.468
    2.019	9.986
    2.088	9.617
    2.192	9.39
    2.314	9.56
    2.366	9.702
    2.418	0.965
    2.816	0.993
    3.18	0.993
    3.544	1.021
    3.891	1.021
    4.289	1.021
    4.671	1.021
    5.009	1.021  ]

    xv      = auxm[:,1]
    pv      = auxm[:,2]

    #Interpolate:
    interp  = LinearInterpolation(xv, pv)
    return interp(xm[:,1])

end 

function ShuOsher_rho_ENO3_N200(xm::Matrix{Float64})

    #Load data:
    auxm,   = readdlm("$(@__DIR__)/ShuOsher_ENO3_N200.csv",',', header=true)
    xv      = auxm[:,1]
    rhov    = auxm[:,2]

    #Interpolate:
    interp  = LinearInterpolation(xv, rhov)
    return interp(xm[:,1])

end

function ShuOsher_rho_ENO3_N400(xm::Matrix{Float64})

    #Load data:
    auxm,   = readdlm("$(@__DIR__)/ShuOsher_ENO3_N400.csv",',', header=true)
    xv      = auxm[:,1]
    rhov    = auxm[:,2]

    #Interpolate:
    interp  = LinearInterpolation(xv, rhov)
    return interp(xm[:,1])

end

function ShuOsher_rho_MUSCL2_N800(xm::Matrix{Float64})

    #Load data:
    auxm,   = readdlm("$(@__DIR__)/ShuOsher_MUSCL2_N800.csv",',', header=true)
    xv      = auxm[:,1]
    rhov    = auxm[:,2]

    #Interpolate:
    interp  = LinearInterpolation(xv, rhov)
    return interp(xm[:,1])

end

function ShuOsher_rho_FV(xm::Matrix{Float64}, N::Int, r::Int)

    @warn "Deprecated"

    #Load data:
    flag    = false
    if r==3
        if N==200
            auxm,   = readdlm("$(@__DIR__)/ShuOsher_ENO3_N200.csv",',', header=true)
            flag    = true
        elseif N==400
            auxm,   = readdlm("$(@__DIR__)/ShuOsher_ENO3_N400.csv",',', header=true)
            flag    = true
        elseif N==800
            @warn "This case is r=2"
            auxm,   = readdlm("$(@__DIR__)/ShuOsher_MUSCL2_N800.csv",',', header=true)
            flag    = true
        end
    else
        auxm    = readdlm("$(@__DIR__)/ShuOsher_N$(N)r$(r).csv", ',')
        flag    = true
    end
    if !flag
        error("Case N=$N, r=$r not defined")
    end

    #Extract and complete data:
    xv      = vcat(-5.05, auxm[:,1], 5.05)
    rhov    = vcat(3.857143, auxm[:,2], 1.0+0.2*sin(5*5.05))

    #Interpolate:
    interp  = LinearInterpolation(xv, rhov)
    return interp(xm[:,1])

end

function ShuOsher_FV(xm::Matrix{Float64}, vble::String, N::Int, r::String)

    #Load data:
    flag    = false
    if r=="3ShuOsher"
        if N==200
            auxm,   = readdlm("$(@__DIR__)/ShuOsher_ENO3_N200.csv",',', header=true)
            flag    = true
        elseif N==400
            auxm,   = readdlm("$(@__DIR__)/ShuOsher_ENO3_N400.csv",',', header=true)
            flag    = true
        end
    elseif r=="5Diaz"
        auxm    = readdlm("$(@__DIR__)/ShuOsher_N$(N)r5.csv", ',')
        flag    = true
    elseif r=="2PhD"
        #MUSCL--Hancock from PhD course
        auxm    = readdlm("$(@__DIR__)/ShuOsher_N$(N)r2.csv", ',')
        flag    = true
    elseif r=="2Diaz"
        #MUSCL from Diaz:
        auxm    = readdlm("$(@__DIR__)/ShuOsher_N$(N)r2bis.csv", ',')
        flag    = true
    elseif r=="2ShuOsher" && N==800
        #From Shu--Osher paper, N=800
        auxm,   = readdlm("$(@__DIR__)/ShuOsher_MUSCL2_N800.csv",',', header=true)
        flag    = true
    end
    if !flag
        error("Case N=$N, r=$r not defined")
    end

    #Select vble:
    if vble=="rho"
        ivar    = 1
        v1      = 3.857143
        v2      = 1.0+0.2*sin(5*5.05)
    elseif vble=="u"
        ivar    = 2
        v1      = 2.629369
        v2      = 0.0
    elseif vble=="p"
        ivar    = 3
        v1      = 10.33333
        v2      = 1.0
    else
        error("Solution for $vble not available")
    end

    #Extract and complete data:
    xv      = vcat(-5.05, auxm[:,1], 5.05)
    vblev   = vcat(v1, auxm[:,1+ivar], v2)

    #Interpolate:
    interp  = LinearInterpolation(xv, vblev)
    return interp(xm[:,1])

end

function FVSol(xm::Matrix{Float64}, vble::String, StudyCase::String, N::Int, r::Int)

    #Load data:
    auxm    = readdlm("$(@__DIR__)/$(StudyCase)_N$(N)r$(r).csv",',')
    
    #Extract and complete data:
    if vble=="rho"
        xv      = auxm[:,1]
        vblev   = auxm[:,2]
    elseif vble=="v"
        xv      = auxm[:,1]
        vblev   = auxm[:,3]
    elseif vble=="p"
        xv      = auxm[:,1]
        vblev   = auxm[:,4]
    else
        error("Solution for $vble not available")
    end

    #Extract and complete data:
    h       = mean(diff(xv))
    xv      = vcat(xv[1]-0.01*h, xv, xv[end]+0.01*h)
    vblev   = vcat(vblev[1], vblev, vblev[end])

    #Interpolate:
    interp  = LinearInterpolation(xv, vblev)
    return interp(xm[:,1])

end

function FV2DSol(StudyCase::String, N::Int, r::Int)

    #Load data:
    auxm    = readdlm("$(@__DIR__)/$(StudyCase)_N$(N)_r$(r).csv",',')
    
    return auxm[:,1], auxm[:,2], auxm[:,3], auxm[:,4] #r, rho, |v|, p

end
