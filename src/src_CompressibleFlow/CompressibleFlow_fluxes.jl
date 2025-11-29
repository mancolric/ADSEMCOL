#-------------------------------------------------------------------------------
#AUXILIARY FUNCTIONS FOR DOMAIN INTEGRALS:

#Update flux terms with contribution of hyperbolic terms. 
#This function receives several matrices: u1,...,uN,p. Then, updates flux_qp,
#which is a nVars x 2 matrix of Matrix{Float64} with the corresponding values of 
#f_{Ij}, I=1,..,,nVars, j=1:2.
#The function also updates dflux_du_qp, which is a nVars x 2 x nVars array with the values of
#df_{Ij}/du_J.
function HyperbolicFlux!(model::GasModel,
    u::Vector{MFloat}, udep::Vector{Vector{MFloat}},
    ComputeJ::Bool,
    flux::Matrix{MFloat},
    dflux::Array{MFloat,3}) where MFloat<:Matrix{Float64}
    
    nSpecies                    = model.nSpecies
    
    #Extract independent variables:
    rho     = udep[DepVarIndex(model,"rho")][1]
    p       = udep[DepVarIndex(model,"p")][1]
    vx      = udep[DepVarIndex(model,"vx")][1]
    vy      = udep[DepVarIndex(model,"vy")][1]
    gamma   = udep[DepVarIndex(model,"gamma")][1]
    
    #Continuity eqs: f_sj = v_j * rhoY_s 
    for ss=1:nSpecies
        @mlv    flux[ss,1]      += vx*u[ss]
        @mlv    flux[ss,2]      += vy*u[ss]
    end
    
    #Momentum eqs: f_{nSpecies+i,j} = v_j*m_i + p*delta_ij
    @mlv    flux[nSpecies+1,1]  += vx*u[nSpecies+1] + p
    @mlv    flux[nSpecies+1,2]  += vy*u[nSpecies+1]
    @mlv    flux[nSpecies+2,1]  += vx*u[nSpecies+2]
    @mlv    flux[nSpecies+2,2]  += vy*u[nSpecies+2] + p
    
    #Energy eqs: f_{nVars,j} = v_j*rhoE + p*v_j
    @mlv    flux[nSpecies+3,1]  += vx*u[nSpecies+3] + p*vx
    @mlv    flux[nSpecies+3,2]  += vy*u[nSpecies+3] + p*vy
        
    #Jacobian terms:
    if ComputeJ
        
        #From 
        #   rhoE    = p/(gamma-1) + 0.5*((rhovx)^2 + (rhovy)^2)/rho,
        #we obtain
        #   0       = 1/(gamma-1) * dp/drho - 0.5*v^2
        #   0       = 1/(gamma-1) * dp/drhovx + vx
        #   0       = 1/(gamma-1) * dp/drhovy + vy
        #   1       = 1/(gamma-1) * dp/drhoE 
        
        #Derivatives of pressure w.r.t. conservative variables:
        dp_du       = udep[DepVarIndex(model,"dp_du")]
        dp_drhoY    = dp_du[1:nSpecies]
        dp_drhovx   = dp_du[nSpecies+1]
        dp_drhovy   = dp_du[nSpecies+2]
        dp_drhoE    = dp_du[nSpecies+3]
        
        #Derivatives of velocities w.r.t. conservative variables:
        #v_i = rhov_i/rho
        dvx_drho    = @mlv -vx/rho 
        dvx_drhovx  = @mlv 1.0/rho
        dvy_drho    = @mlv -vy/rho 
        dvy_drhovy  = @mlv 1.0/rho
        
        #Continuity eqs:
        for s1=1:nSpecies
            for s2=1:nSpecies 
                @mlv    dflux[s1,1,s2]  += dvx_drho*u[s1] + vx*(s1==s2)
                @mlv    dflux[s1,2,s2]  += dvy_drho*u[s1] + vy*(s1==s2)
            end
            @mlv    dflux[s1,1,nSpecies+1]  += dvx_drhovx*u[s1]
            @mlv    dflux[s1,2,nSpecies+2]  += dvy_drhovy*u[s1]
        end
        
        #x1-momentum equation: 
        II                              = nSpecies+1
        for ss=1:nSpecies
            @mlv    dflux[II,1,ss]      += dvx_drho*u[nSpecies+1] + dp_drhoY[ss]
            @mlv    dflux[II,2,ss]      += dvy_drho*u[nSpecies+1]
        end
        @mlv    dflux[II,1,nSpecies+1]  += dvx_drhovx*u[nSpecies+1] + vx + dp_drhovx
        @mlv    dflux[II,1,nSpecies+2]  += dp_drhovy
        @mlv    dflux[II,1,nSpecies+3]  += dp_drhoE
        @mlv    dflux[II,2,nSpecies+1]  += vy
        @mlv    dflux[II,2,nSpecies+2]  += dvy_drhovy*u[nSpecies+1]
        
        #x2-momentum equation: 
        II                              = nSpecies+2
        for ss=1:nSpecies
            @mlv    dflux[II,1,ss]      += dvx_drho*u[nSpecies+2]
            @mlv    dflux[II,2,ss]      += dvy_drho*u[nSpecies+2] + dp_drhoY[ss]
        end
        @mlv    dflux[II,1,nSpecies+1]  += dvx_drhovx*u[nSpecies+2]
        @mlv    dflux[II,1,nSpecies+2]  += vx
        @mlv    dflux[II,2,nSpecies+1]  += dp_drhovx
        @mlv    dflux[II,2,nSpecies+2]  += dvy_drhovy*u[nSpecies+2] + vy + dp_drhovy
        @mlv    dflux[II,2,nSpecies+3]  += dp_drhoE
        
        #Energy equation: f_I1=p*rhovx/rho, f_I2=p*rhovy/rho
        II                              = nSpecies+3
        for ss=1:nSpecies
            @mlv    dflux[II,1,ss]      += dvx_drho*u[nSpecies+3] + dp_drhoY[ss]*vx - p*vx/rho
            @mlv    dflux[II,2,ss]      += dvy_drho*u[nSpecies+3] + dp_drhoY[ss]*vy - p*vy/rho
        end
        @mlv    dflux[II,1,nSpecies+1]  += dvx_drhovx*u[nSpecies+3] + dp_drhovx*vx + p/rho
        @mlv    dflux[II,1,nSpecies+2]  += dp_drhovy*vx
        @mlv    dflux[II,1,nSpecies+3]  += vx + dp_drhoE*vx
        @mlv    dflux[II,2,nSpecies+1]  += dp_drhovx*vy
        @mlv    dflux[II,2,nSpecies+2]  += dvy_drhovy*u[nSpecies+3] + dp_drhovy*vy + p/rho
        @mlv    dflux[II,2,nSpecies+3]  += vy + dp_drhoE*vy
        
    end

    return
    
end

#Update flux terms with contribution of viscous flux:
function ViscousFlux!(model::GasModel, 
    nu::MFloat, beta::MFloat,
    u::Vector{MFloat}, udep::Vector{Vector{MFloat}}, du::Matrix{MFloat}, 
    ComputeJ::Bool, 
    flux::Matrix{MFloat}, dflux_du::Array{MFloat,3},  dflux_dgradu::Array{MFloat,4}) where MFloat<:Matrix{Float64}
        
    nSpecies        = model.nSpecies
    
    #Extract variables:
    rhov            = [u[nSpecies+2],u[nSpecies+3]]
    rho             = udep[DepVarIndex(model,"rho")][1] #MFloat
    vx              = udep[DepVarIndex(model,"vx")][1]
    vy              = udep[DepVarIndex(model,"vy")][1]
    v               = [vx, vy]
    
    #Extract gradients of rho:
    drho            = [ sum(du[1:nSpecies,1]), sum(du[1:nSpecies,2]) ]
    
    #Derivatives of some variables:
    dv_duJ          = Matrix{MFloat}(undef,2,3)     #Derivs w.r.t. rho, rhovx, rhovy:
    if ComputeJ
        for ii=1:2
            #v_i = rhov_i/rho;
            dv_duJ[ii,1]    = @mlv    -v[ii]/rho
            dv_duJ[ii,2]    = @mlv    (ii==1)/rho
            dv_duJ[ii,3]    = @mlv    (ii==2)/rho
        end
    end
    
    #---------------------------------------------------------------------
    #Action of viscous flux:
    
    #f_{nSpecies+i,j}   = -sigma_ij
    #f_{nSpecies+3,j}   = -v_i sigma_ij
    #sigma_ij           = rho*nu*(d_j v_i + d_i v_j) + rho*(beta-2/3*nu) d_k v_k delta_ij
    #rho d_i v_j        = d_i(rho v_j) - v_j d_i rho
    
    #Divergence term:
    rho_divv        = @mlv    du[nSpecies+1,1]-v[1]*drho[1] + du[nSpecies+2,2]-v[2]*drho[2]
    drho_divv_duJ   = Vector{MFloat}(undef,3)   #Derivatives w.r.t. rho, rhovx, rhovy
    if ComputeJ
        for JJ=1:3
            drho_divv_duJ[JJ]       = @mlv    -dv_duJ[1,JJ]*drho[1]-dv_duJ[2,JJ]*drho[2]
        end
    end
    
    #Allocate sigma and loop:
    sigma_ij                        = zeros(size(flux[1,1]))
    for jj=1:2, ii=1:2
    
        @mlv    sigma_ij            = nu*(  du[nSpecies+ii,jj]-v[ii]*drho[jj] + 
                                            du[nSpecies+jj,ii]-v[jj]*drho[ii] ) + 
                                      (beta-2/3*nu) * rho_divv * (ii==jj)
        @mlv    flux[nSpecies+ii,jj]-= sigma_ij
        @mlv    flux[nSpecies+3,jj] -= v[ii]*sigma_ij
        
        if ComputeJ
        
            #Derivatives w.r.t. the gradients:
            for ll=1:2
                for ss=1:nSpecies #loop species
                    @mlv    dflux_dgradu[nSpecies+ii,jj,ss,ll]  -=
                        nu*( -v[ii]*(jj==ll) - v[jj]*(ii==ll) ) +
                        (beta-2/3*nu) * (-v[1]*(ll==1)-v[2]*(ll==2)) * (ii==jj)
                    @mlv    dflux_dgradu[nSpecies+3,jj,ss,ll]   -= 
                        v[ii] * (
                            nu*( -v[ii]*(jj==ll) - v[jj]*(ii==ll) ) +
                            (beta-2/3*nu) * (-v[1]*(ll==1)-v[2]*(ll==2)) * (ii==jj) )
                end
                for mm=1:2 #loop momentum
                    @mlv    dflux_dgradu[nSpecies+ii,jj,nSpecies+mm,ll] -= 
                        nu*( (ii==mm)*(jj==ll) + (jj==mm)*(ii==ll) ) + 
                        (beta-2/3*nu)* ((mm==1)*(ll==1) + (mm==2)*(ll==2)) * (ii==jj)
                    @mlv    dflux_dgradu[nSpecies+3,jj,nSpecies+mm,ll]  -= 
                        v[ii] * (
                            nu*( (ii==mm)*(jj==ll) + (jj==mm)*(ii==ll) ) + 
                            (beta-2/3*nu)* ((mm==1)*(ll==1) + (mm==2)*(ll==2)) * (ii==jj) )
                end
            end
        
            #Derivatives w.r.t. conservative variables are ignored:
            #=
            for ss=1:nSpecies #loop species
                @mlv    dflux_du[nSpecies+ii,jj,ss]             -=
                    nu*(    -dv_duJ[ii,1]*drho[jj] + 
                            -dv_duJ[jj,1]*drho[ii] ) + 
                    (beta-2/3*nu) * drho_divv_duJ[1] * (ii==jj)
                @mlv    dflux_du[nSpecies+3,jj,ss]              -=
                    v[ii] *( nu*(   -dv_duJ[ii,1]*drho[jj] + 
                                    -dv_duJ[jj,1]*drho[ii] ) + 
                            (beta-2/3*nu) * drho_divv_duJ[1] * (ii==jj) ) + 
                    dv_duJ[ii,1] * sigma_ij
            end
            for mm=1:2  #loop momentum
                @mlv    dflux_du[nSpecies+ii,jj,nSpecies+mm]    -= 
                    nu*(    -dv_duJ[ii,1+mm]*drho[jj] +
                            -dv_duJ[jj,1+mm]*drho[ii] ) +
                    (beta-2/3*nu) * drho_divv_duJ[1+mm] * (ii==jj) 
                @mlv    dflux_du[nSpecies+3,jj,nSpecies+mm]     -= 
                    v[ii] *(nu*(    -dv_duJ[ii,1+mm]*drho[jj] +
                                    -dv_duJ[jj,1+mm]*drho[ii] ) +
                            (beta-2/3*nu) * drho_divv_duJ[1+mm] * (ii==jj) ) + 
                    dv_duJ[ii,1+mm] * sigma_ij
            end
            =#
            
        end
    
    end
    
    return
    
end
    
function HeatFlux!(model::GasModel, 
    kappa_rho_cv::MFloat,
    u::Vector{MFloat}, udep::Vector{Vector{MFloat}}, du::Matrix{MFloat}, 
    ComputeJ::Bool, 
    flux::Matrix{MFloat}, dflux_du::Array{MFloat,3},  dflux_dgradu::Array{MFloat,4}) where MFloat<:Matrix{Float64}
        
    nSpecies        = model.nSpecies
    
    #Extract variables:
    rhov            = [u[nSpecies+2],u[nSpecies+3]]
    rho             = udep[DepVarIndex(model,"rho")][1] #MFloat
    vx              = udep[DepVarIndex(model,"vx")][1]
    vy              = udep[DepVarIndex(model,"vy")][1]
    v               = [vx, vy]
    e_i             = udep[DepVarIndex(model,"e_i")]    #Vector{MFloat}
    
    #Extract gradients of rho:
    drho            = [ sum(du[1:nSpecies,1]), sum(du[1:nSpecies,2]) ]
    
    #---------------------------------------------------------------------
    #Action of heat transfer. Note that de = cv*dT + sum_k e_k dY_k
    
    # q_j               = -kappa d_j T 
    # rho e             = rho Y_k e_k 
    # d(rho e)          = d(rho Y_k) e_k + rho Y_k de_k 
    #                   = d(rho Y_k) e_k + rho Y_k cv_k dT
    #                   = d(rho Y_k) e_k + rho cv dT
    # dT                = 1/rho/cv * (d(rho e) - e_k d(rho Y_k))
    # rho e             = rho E - 0.5 * ( (rhovx)^2 + (rhovy)^2 ) / rho
    # d(rho e)          = d(rho E) - rhovx/rho * d(rhovx) - rhovy/rho * d(rhovy) + 0.5*((rhovx)^2 + (rhovy)^2)/rho^2 * d(rho) 
    for jj=1:2
    
        @mlv    flux[nSpecies+3,jj]     -= kappa_rho_cv * 
            (du[nSpecies+3,jj] - (vx*du[nSpecies+1,jj] + vy*du[nSpecies+2,jj]) + 
                0.5*(vx*vx+vy*vy)*drho[jj])
        for ss=1:nSpecies
            @mlv    flux[nSpecies+3,jj] += kappa_rho_cv * (e_i[ss] * du[ss,jj])
        end
        
        if ComputeJ
            
            #Derivatives w.r.t. gradients:
            for ss=1:nSpecies
                @mlv    dflux_dgradu[nSpecies+3,jj,ss,jj]       -= 
                    kappa_rho_cv*(0.5*(vx*vx+vy*vy) - e_i[ss])
            end
            @mlv    dflux_dgradu[nSpecies+3,jj,nSpecies+1,jj]   -= 
                kappa_rho_cv * (-vx)
            @mlv    dflux_dgradu[nSpecies+3,jj,nSpecies+2,jj]   -= 
                kappa_rho_cv * (-vy)
            @mlv    dflux_dgradu[nSpecies+3,jj,nSpecies+3,jj]   -= 
                kappa_rho_cv
        
            #Derivatives w.r.t. variables are ignored.
            #=
            for ss=1:nSpecies
                @mlv    dflux_du[nSpecies+3,jj,ss]  -= kappa_rho_cv * (
                    -dv[1]*du[nSpecies+1,jj]-dv[2]*du[nSpecies+2,jj] + 
                    (vx*dv[1] + vy*dv[2])*drho[jj] )
            end
            =#
            
        end
            
    end
    
    return
    
end

function BrennerMassFlux!(model::GasModel, 
    epsilon::MFloat,
    u::Vector{MFloat}, udep::Vector{Vector{MFloat}}, du::Matrix{MFloat}, 
    ComputeJ::Bool, 
    flux::Matrix{MFloat}, dflux_du::Array{MFloat,3},  dflux_dgradu::Array{MFloat,4}) where MFloat<:Matrix{Float64}
        
    nSpecies        = model.nSpecies
    
    #Each specie k has a diffusion velocity -epsilon d_j (rhoY_k)/rhoY_k.
    #That is, the flux of the variable phi (per mass unit) associated 
    #with mass diffusion of k specie is
    #   f_j     = v_j rhoY_k phi = -epsilon d_j(rhoY_k) phi, 
    #with phi=1, v_i, E_i.
    #That is, for the mass equations,
    #   f_j     = -epsilon d_j(rhoY_k).
    #For the momentum equations,
    #   f_j     = sum_k [-epsilon d_j(rhoY_k) v_i] 
    #           = -epsilon (d_j rho) v_i
    #For the total energy equation,
    #   f_j     = sum_k (-epsilon d_j(rhoY_k))(e_i + 0.5 v^2)
    #           = sum_k (-epsilon d_j(rhoY_k))e_i + 
    #               -epsilon (d_j rho) * 0.5 v^2 
    
    #Extract variables:
    rho             = udep[DepVarIndex(model,"rho")][1] #MFloat
    vx              = udep[DepVarIndex(model,"vx")][1]
    vy              = udep[DepVarIndex(model,"vy")][1]
    v               = [vx, vy]
    e_i             = udep[DepVarIndex(model,"e_i")]    #Vector{MFloat}
    
    #Extract gradients of rho:
    drho            = [ sum(du[1:nSpecies,1]), sum(du[1:nSpecies,2]) ]
    
    for jj=1:2
        
        for ss=1:nSpecies
            @mlv    flux[ss,jj]         -= epsilon*du[ss,jj]
        end
        for II=nSpecies+1:nSpecies+2
            @mlv    flux[II,jj]         -= epsilon*drho[jj]*u[II]/rho
        end
        for II=nSpecies+3
            @mlv flux[II,jj]    -= epsilon*drho[jj]*0.5*(vx*vx+vy*vy)
            for ss=1:nSpecies
                @mlv flux[II,jj]    -= epsilon*du[ss,jj]*e_i[ss]
            end
        end
        
        if ComputeJ
        
            #Derivatives w.r.t. gradients of conservative variables:
            for ss=1:nSpecies
                @mlv    dflux_dgradu[ss,jj,ss,jj]       -= epsilon
                for II=nSpecies+1:nSpecies+2
                    @mlv    dflux_dgradu[II,jj,ss,jj]   -= epsilon*u[II]/rho
                end
                for II=nSpecies+3
                    @mlv dflux_dgradu[II,jj,ss,jj]      -= epsilon*(e_i[ss]+0.5*(vx*vx+vy*vy)) 
                end
            end
            
            #Derivatives w.r.t. conservative variables are ignored
            
        end
        
    end
    
    return
    
end

function MassDiffusionFlux!(model::GasIdeal, 
    u::Vector{MFloat}, udep::Vector{Vector{MFloat}}, du::Matrix{MFloat}, 
    ComputeJ::Bool, 
    flux::Matrix{MFloat}, dflux_du::Array{MFloat,3},  dflux_dgradu::Array{MFloat,4}) where MFloat<:Matrix{Float64}

    return
    
end

function MassDiffusionFlux!(model::ReactiveGas, 
    u::Vector{MFloat}, udep::Vector{Vector{MFloat}}, du::Matrix{MFloat}, 
    ComputeJ::Bool, 
    flux::Matrix{MFloat}, dflux_du::Array{MFloat,3},  dflux_dgradu::Array{MFloat,4}) where MFloat<:Matrix{Float64}
        
    nSpecies        = model.nSpecies
    
    #Each specie k has a diffusion velocity (v_k)_j = - D_k d_j(Y_k) / Y_k.
    #with 
    #   sum_k [(v_k)_j rhoY_k] = 0
    #Note that if all the species have the same coefficient D, then the latter 
    #condition is automatically satisfied.
    #
    #Hence, for each mass equation we have to add a flux
    #   f_j     = (v_k)_j rhoY_k = -D_k d_j(Y_k) * rho
    #For the momentum equations, we add
    #   f_j     = sum_k [(v_k)_j rhoY_k v_i] = 0
    #For the energy equations, we add
    #   f_j     = sum_k [(v_k)_j rhoY_k (e_k + 0.5 v^2)] = 
    #           = sum_k [(v_k)_j rhoY_k e_k] + 0
    #which is the generalized Fourier law.
    #
    #Finally, note that
    #   J_kj    = (v_k)_j rhoY_k = -D_k d_j(Y_k) * rho 
    #           = -D_k [d_j(rhoY_k) - Y_k d_j(rho)]
    #is the so-called mass flux, and then
    #   f_mass_k_j  = J_kj
    #   f_energy_j  = sum_k J_kj * e_k
    
    #Extract variables:
    rho             = udep[DepVarIndex(model,"rho")][1] #MFloat
    vx              = udep[DepVarIndex(model,"vx")][1]
    vy              = udep[DepVarIndex(model,"vy")][1]
    v               = [vx, vy]
    h_i             = udep[DepVarIndex(model,"h_i")]    #Vector{MFloat}
    e_i             = udep[DepVarIndex(model,"e_i")]    #Vector{MFloat}
    D_i             = udep[DepVarIndex(model,"D_i")]    #Vector{MFloat}

    #Extract gradients of rho:
    drho            = [ sum(du[1:nSpecies,1]), sum(du[1:nSpecies,2]) ]
    
    for jj=1:2
        
        for s1=1:nSpecies
            @mlv flux[s1,jj]            -= D_i[s1]*(du[s1,jj]-u[s1]/rho*drho[jj])
            @mlv flux[nSpecies+3,jj]    -= D_i[s1]*(du[s1,jj]-u[s1]/rho*drho[jj]) * 
                                            h_i[s1]
        end
        
        if ComputeJ
        
            #Derivatives w.r.t. gradients of conservative variables:
            for s1=1:nSpecies, s2=1:nSpecies
                @mlv dflux_dgradu[s1,jj,s2,jj]          -= D_i[s1]*((s1==s2)-u[s1]/rho)
                @mlv dflux_dgradu[nSpecies+3,jj,s2,jj]  -= D_i[s1]*((s1==s2)-u[s1]/rho)*
                                                            h_i[s1]
            end
            
            #Derivatives w.r.t. conservative variables are ignored
            
        end
        
    end
    
    return
    
end

#Add monolithic diffusion:
#   f_Ij = -tau d_j u_I
function epsilonFlux!(model::GasModel, tau::MFloat, duB::Matrix{MFloat}, 
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

function TestFlux!(model::GasIdeal, 
    epsilon::MFloat, nu::MFloat, beta::MFloat, kappa_rho_cv::MFloat,
    u::Vector{MFloat}, udep::Vector{Vector{MFloat}}, du::Matrix{MFloat}, 
    ComputeJ::Bool, 
    flux::Matrix{MFloat}, dflux_du::Array{MFloat,3}, dflux_dgradu::Array{MFloat,4}) where MFloat<:Matrix{Float64}
        
    nSpecies            = model.nSpecies
    
    for II=1:1, jj=1:2
        @mlv    flux[II,jj]                     -= nu*du[II,jj]
        if ComputeJ
            @mlv    dflux_dgradu[II,jj,II,jj]   -= nu
        end
    end
    
    for II=1:nSpecies+3, jj=1:2
        @mlv    flux[II,jj]                     -= beta*u[II]
        if ComputeJ
            @mlv    dflux_du[II,jj,II]          -= beta
        end
    end
    
    return
    
end

function TestFlux!(model::GasIdeal, 
    u::Vector{MFloat}, udep::Vector{Vector{MFloat}}, ComputeJ::Bool, 
    flux::Matrix{MFloat}, dflux_du::Array{MFloat,3}) where MFloat<:Matrix{Float64}
        
    nSpecies            = model.nSpecies
    
    beta                = 2.0
    for II=1:nSpecies+3, jj=2
        @mlv    flux[II,jj]                     += beta*u[II]
        if ComputeJ
            @mlv    dflux_du[II,jj,II]          += beta
        end
    end
    
    return
    
end

#Function that evaluates the source term:
function source!(model::GasIdeal, t::Float64, x::Vector{MFloat}, 
    u::Vector{MFloat}, udep::Vector{Vector{MFloat}}, 
    ComputeJ::Bool, flux::Vector{MFloat}, dflux_du::Matrix{MFloat}) where MFloat<:Matrix{Float64}
    
    nSpecies        = model.nSpecies
    g               = model.g
    if g==0.0
        return
    end
    
    #Extract variables:
    rho             = udep[DepVarIndex(model,"rho")][1] #MFloat
    rhovy           = u[nSpecies+2]
    
    #Gravity terms:
    @avxt @. flux[nSpecies+2]   -= rho*g
    @avxt @. flux[nSpecies+3]   -= rhovy*g
    if ComputeJ
        for ss=1:nSpecies
            @avxt @. dflux_du[nSpecies+2, ss]       -= g
        end
        @avxt @. dflux_du[nSpecies+3, nSpecies+2]   -= g
    end
    
    return
    
end

function source!(model::ReactiveGas, t::Float64, x::Vector{MFloat}, 
    u::Vector{MFloat}, udep::Vector{Vector{MFloat}}, 
    ComputeJ::Bool, flux::Vector{MFloat}, dflux_du::Matrix{MFloat}) where MFloat<:Matrix{Float64}
    
    nSpecies    = model.nSpecies
    mdot_i      = udep[DepVarIndex(model, "mdot_i")] #Vector{MFloat}
    for ss=1:nSpecies
        @mlv flux[ss]               += mdot_i[ss]
    end
    if ComputeJ
        dmdot_ij    = reshape(udep[DepVarIndex(model, "dmdot_ij")], nSpecies, nSpecies+3)   #Matrix{Float64}
        for ss=1:nSpecies, JJ=1:nSpecies+3
            @mlv dflux_du[ss,JJ]    += dmdot_ij[ss,JJ]
        end
    end
    
    return
    
end

