#--------------------------------------------------------------------------------
#BOUNDARY CONDITIONS:

#Impose penalty flux fn_I = sigma*(u_I - uBC_I):
function penalty!(model::SCD, 
    sigma::MFloat,
    u::Vector{MFloat}, uBC::Vector{MFloat}, duBC_du::Matrix{MFloat}, 
    ComputeJ::Bool,
    bflux::Vector{MFloat}, dbflux_du::Matrix{MFloat}) where MFloat<:Matrix{Float64}
    
    #Update boundary flux:
    for II=1:model.nVars
        @mlv    bflux[II]                   += sigma*(u[II]-uBC[II])
        if ComputeJ
            for JJ=1:model.nVars
                @mlv    dbflux_du[II,JJ]    += sigma*((II==JJ)-duBC_du[II,JJ])
            end
        end
    end
    
    return
    
end

#Dirichlet conditions:
function bflux!(model::SCD, BC::Dirichlet, _bqp::TrBintVars, ComputeJ::Bool)
    
    nb                      = _bqp.nb
    t                       = _bqp.t
    x                       = _bqp.x
    u                       = _bqp.u
    du                      = _bqp.gradu
    ParentElems             = _bqp.Binteg2D.bmesh.ParentElems
    metric                  = _bqp.Binteg2D.mesh.metric
    
    #Dirichlet condition:
    uBC                     = BC.uDir(t,x)
    duBC_du                 = Matrix{Matrix{Float64}}(undef,1,1)
    duBC_du[1,1]            = zeros(size(u[1]))
    
    #=
    #Hyperbolic term:
    a                               = model.a(t,x,u)
    an                              = @mlv a[1]*nb[1]+a[2]*nb[2]
    inflow                          = @mlv an <= 0.0
    outflow                         = @mlv !inflow
    @mlv _bqp.f[1][inflow]          += an[inflow]*uBC[1][inflow]
    @mlv _bqp.f[1][outflow]         += an[outflow]*u[1][outflow]
    if ComputeJ
        da_du                       = model.da_du(t,x,u)
        dan_du                      = @. da_du[1]*nb[1] + da_du[2]*nb[2]
        @. _bqp.df_du[1,1][inflow]  += dan_du[inflow]*uBC[1][inflow]
        @. _bqp.df_du[1,1][outflow] += dan_du[outflow]*u[1][outflow] + an[outflow]
    end
    =#
    
    #Hyperbolic term:
    a                           = model.a(t,x,u)
    a_n                         = @avxt @. a[1]*nb[1]+a[2]*nb[2]
    da_du                       = model.da_du(t,x,u)
    ahat_1                      = @avxt @. da_du[1]*u[1] + a[1]
    ahat_2                      = @avxt @. da_du[2]*u[1] + a[2]
    ahat_n                      = @avxt @. ahat_1*nb[1] + ahat_2*nb[2]
    inflow                      = @avxt @. ahat_n <= 0.0
    outflow                     = @avxt @. !inflow
    #
    aBC                         = model.a(t,x,uBC)
    aBC_n                       = @avxt @. aBC[1]*nb[1] + aBC[2]*nb[2]
    @avxt @. _bqp.f[1][inflow]  += aBC_n[inflow]*uBC[1][inflow]
    @avxt @. _bqp.f[1][outflow] += a_n[outflow]*u[1][outflow]
    if ComputeJ
        @avxt @. _bqp.df_du[1,1][inflow]    += 0.0
        @avxt @. _bqp.df_du[1,1][outflow]   += ahat_n[outflow]
    end
    
    #Allocate viscous fluxes:
    flux, dflux_du, dflux_dgradu        = FluxAllocate(1, size(u[1]), ComputeJ)
    epsilon, depsilon_du, 
        depsilon_dgradu                 = ViscosityAllocate(1, size(u[1]), ComputeJ)
    
    #Natural viscosity:
    epsilon                             = model.epsilon(t,x,u)[1]
    depsilon_du                         = zeros(0,0)
    if ComputeJ
        depsilon_du                     = model.depsilon_du(t,x,u)[1]
    end
    
    #Extrapolate natural viscous flux:
    DiffusiveFlux!(model, epsilon, depsilon_du, u, du, ComputeJ, 
        flux, dflux_du, dflux_dgradu)
    fn, dfn_du, dfn_dgradu              = ProjectFlux(flux, dflux_du, dflux_dgradu, nb, ComputeJ)     
    #
    @mlv _bqp.f[1]                  += fn[1]
    if ComputeJ
        @mlv _bqp.df_du[1,1]        += dfn_du[1,1]
        @mlv _bqp.df_dgradu[1,1,1]  += dfn_dgradu[1,1,1]
        @mlv _bqp.df_dgradu[1,1,2]  += dfn_dgradu[1,1,2]
    end
    #Add penalty terms:
    h_Elems                 = @mlv 1.0/sqrt(metric.lambda_bar[ParentElems])
    hp                      = h_Elems./_bqp.FesOrder * ones(1, _bqp.nqp)
    sigma                   = @mlv model.CW * epsilon / hp
    penalty!(model, sigma, u, uBC, duBC_du, ComputeJ, _bqp.f, _bqp.df_du)
    
    return
    
end

#Neumann conditions:
function bflux!(model::SCD, BC::Neumann, _bqp::TrBintVars, ComputeJ::Bool)
    
    nb                      = _bqp.nb
    t                       = _bqp.t
    x                       = _bqp.x
    u                       = _bqp.u

    #Hyperbolic term:
    a                       = model.a(t,x,u)
    @mlv _bqp.f[1]          += (a[1]*nb[1]+a[2]*nb[2])*u[1]
    if ComputeJ
        da_du               = model.da_du(t,x,u)
        @mlv _bqp.df_du[1,1]+= (da_du[1]*nb[1]+da_du[2]*nb[2])*u[1] + 
                                a[1]*nb[1]+a[2]*nb[2]
    end
    
    #Since u is extrapolated, penalty is zero.
    
    #Impose viscous flux:
    fn                      = BC.q(t,x)
    @mlv _bqp.f[1]          += fn[1]
    
    return
    
end
