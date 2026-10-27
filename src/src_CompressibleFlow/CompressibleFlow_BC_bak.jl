#-------------------------------------------------------------------------------
#ROTATION MATRICES:

function uRotate(nb::Vector{MFloat}, tb::Vector{MFloat}, u::Vector{MFloat}) where MFloat<:AbstractMatrix{Float64}

    nVars           = length(u)
    udims           = size(u)
    q               = Vector{MFloat}(undef,nVars)
    for II=1:nVars
        q[II]       = zeros(udims)
    end
    uRotate!(nb, tb, u, q)
    return q
    
end
    
function uRotate!(nb::Vector{MFloat}, tb::Vector{MFloat}, u::Vector{MFloat}, q::Vector{MFloat}) where MFloat<:AbstractMatrix{Float64}

    nVars           = length(u)

    #Rotate variables at the boundary. We have u[II]=q[II], except for 
    #   q[nSpecies+1], which is the normal momentum
    #   q[nSpecies+2], which is the tangent momentum
    for II in vcat(1:nVars-3, nVars)
        q[II]       = @mlv    u[II]
    end
    #v_n    = v_x * n_x + v_y * n_y
    q[nVars-2]      = @mlv    u[nVars-2]*nb[1] + u[nVars-1]*nb[2]
    #t_n    = v_x * t_x + v_y * t_y
    q[nVars-1]      = @mlv    u[nVars-2]*tb[1] + u[nVars-1]*tb[2]
    
    return q
    
end

#Matrix to rotate u and obtain q:
function du_dq_boundary(nb::Vector{MFloat}, tb::Vector{MFloat}, nVars::Int) where MFloat<:AbstractMatrix{Float64}

    R                   = Matrix{MFloat}(undef,nVars,nVars)
    udims               = size(nb[1])
    for JJ=1:nVars, II=1:nVars
        R[II,JJ]                = zeros(udims)
    end
    for II in vcat(1:nVars-3, nVars)
        @mlv    R[II,II]        = 1.0
    end
    @mlv    R[nVars-2,nVars-2]  = nb[1]
    @mlv    R[nVars-2,nVars-1]  = nb[2]
    @mlv    R[nVars-1,nVars-2]  = tb[1]
    @mlv    R[nVars-1,nVars-1]  = tb[2]
    
    return R
    
end

function qRotate(nb::Vector{MFloat}, tb::Vector{MFloat}, q::Vector{MFloat}) where MFloat<:AbstractMatrix{Float64}

    nVars           = length(q)
    udims           = size(q)
    u               = Vector{MFloat}(undef,nVars)
    for II=1:nVars
        u[II]       = zeros(udims)
    end
    qRotate!(nb, tb, q, u)
    return u
    
end

function qRotate!(nb::Vector{MFloat}, tb::Vector{MFloat}, q::Vector{MFloat}, u::Vector{MFloat}) where MFloat<:AbstractMatrix{Float64}

    nVars           = length(u)

    #Rotate variables at the boundary. We have u[II]=q[II], except for 
    #   q[nSpecies+1], which is the normal momentum
    #   q[nSpecies+2], which is the tangent momentum:
    for II in vcat(1:nVars-3, nVars)
        u[II]       = @mlv    q[II]
    end
    #v_i    = v_n*n_i + v_t*t_i:
    u[nVars-2]      = @mlv    q[nVars-2]*nb[1] + q[nVars-1]*tb[1]
    u[nVars-1]      = @mlv    q[nVars-2]*nb[2] + q[nVars-1]*tb[2]
    
    return q
    
end

#Matrix to rotate u and obtain q:
function dq_du_boundary(nb::Vector{MFloat}, tb::Vector{MFloat}, nVars::Int) where MFloat<:AbstractMatrix{Float64}

    R                   = Matrix{MFloat}(undef,nVars,nVars)
    udims               = size(nb[1])
    for JJ=1:nVars, II=1:nVars
        R[II,JJ]        = zeros(udims)
    end
    for II in vcat(1:nVars-3, nVars)
        @mlv    R[II,II]        = 1.0
    end
    @mlv    R[nVars-2,nVars-2]  = nb[1]
    @mlv    R[nVars-2,nVars-1]  = tb[1]
    @mlv    R[nVars-1,nVars-2]  = nb[2]
    @mlv    R[nVars-1,nVars-1]  = tb[2]
    
    return R
    
end

#--------------------------------------------------------------------------------
#BOUNDARY CONDITIONS:

#Impose penalty flux fn_I = sigma*(u_I - uBC_I):
function penalty!(model::GasModel, 
    h::MFloat,
    u::Vector{MFloat}, uBC::Vector{MFloat}, duBC_du::Matrix{MFloat}, 
    udep::Vector{Vector{MFloat}}, ComputeJ::Bool,
    bflux::Vector{MFloat}, dbflux_du::Matrix{MFloat}) where MFloat<:Matrix{Float64}
    
    #Viscosities:
    D_penalty       = udep[DepVarIndex(model,"D_penalty")][1]
    sigma           = @mlv model.CW * D_penalty / h
    
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
    
#SlipAdiabatic boundary. 3+nSpecies conditions:
#   m_n             = 0     zero normal velocity
#   fmass_(k,j) n_j = 0     no mass flux, k=1,...,nSpecies
#   tau_tj n_j      = 0     no viscous force in tangent direction
#   q_j n_j         = 0     no heat flux
function qBC_old(model::GasModel, BC::SlipAdiabatic, 
    t::Float64, x::Vector{MFloat}, q::Vector{MFloat},
    ComputeJ::Bool) where MFloat<:Matrix{Float64}
    
    nSpecies                = model.nSpecies
    nVars                   = model.nVars
    
    #Compute qBC and derivatives:
    qBC                     = Vector{MFloat}(undef,nVars)
    for II in vcat(1:nSpecies, nSpecies+2, nSpecies+3)
        qBC[II]             = @mlv    q[II]
    end
    qBC[nSpecies+1]         = @mlv    0.0*q[nSpecies+1]
    dqBC_dq                 = Matrix{MFloat}(undef,nVars,nVars)
    if ComputeJ
        alloc!(dqBC_dq, size(q[1]))
        for II in vcat(1:nSpecies, nSpecies+2, nSpecies+3)
            @mlv    dqBC_dq[II,II]              = 1.0
        end
        @mlv    dqBC_dq[nSpecies+1, nSpecies+1] = 0.0
    end
    
    return qBC, dqBC_dq
    
end
function qBC_(model::GasModel, BC::SlipAdiabatic, 
    t::Float64, x::Vector{MFloat}, q::Vector{MFloat},
    ComputeJ::Bool) where MFloat<:Matrix{Float64}
    
    nSpecies                = model.nSpecies
    nVars                   = model.nVars
    
    #Allocate:
    qBC                     = Vector{MFloat}(undef,nVars)
    dqBC_dq                 = Matrix{MFloat}(undef,nVars,nVars)
    
    #Same density and momentum in tangent direction:
    for II in vcat(1:nSpecies, nSpecies+2)
        qBC[II]             = @tturbo @. q[II]
    end
    if ComputeJ
        alloc!(dqBC_dq, size(q[1]))
        for II in vcat(1:nSpecies, nSpecies+2)
            @tturbo @. dqBC_dq[II,II]   = 1.0
        end
    end
       
    #Zero momentum in normal direction:
    qBC[nSpecies+1]         = @tturbo @. 0.0*q[nSpecies+1]
    
    #Extrapolate internal energy. Note that
    #   rhoE    = rhoe + 0.5*(m_n^2 + m_t^2)/rho
    #   rhoE_BC = rhoe + 0.5*m_t^2/rho
    #so
    #   rhoE_BC = rhoE - 0.5*m_n^2/rho
    xout                = DepVars(model, t, x, q, ["rho"])
    rho                 = xout[1][1]
    qBC[nSpecies+3]     = @tturbo @. q[nSpecies+3] - 0.5*q[nSpecies+1]^2/rho
    if ComputeJ
        for JJ=1:nSpecies
            @tturbo @. dqBC_dq[nSpecies+3, JJ]      = +0.5*q[nSpecies+1]^2/rho^2
        end
        @tturbo @. dqBC_dq[nSpecies+3, nSpecies+1]  = -q[nSpecies+1]/rho
        @tturbo @. dqBC_dq[nSpecies+3, nSpecies+3]  = 1.0
    end
    
    return qBC, dqBC_dq
    
end

#SlipAdiabatic boundary. 3+nSpecies conditions:
#   m_n             = 0     zero normal velocity
#   fmass_(k,j) n_j = 0     no mass flux, k=1,...,nSpecies
#   tau_tj n_j      = 0     no viscous force in tangent direction
#   q_j n_j         = 0     no heat flux
function bflux!(model::GasModel, BC::SlipAdiabatic, 
    _bqp::TrBintVars, ComputeJ::Bool)
    
    nSpecies                = model.nSpecies
    nVars                   = nSpecies+3
    
    nb                      = _bqp.nb
    tb                      = _bqp.tb
    t                       = _bqp.t
    x                       = _bqp.x
    u                       = _bqp.u
    du                      = _bqp.gradu
    ParentElems             = _bqp.Binteg2D.bmesh.ParentElems
    metric                  = _bqp.Binteg2D.mesh.metric
    bflux                   = _bqp.f
    dbflux_du               = _bqp.df_du
    dbflux_dgradu           = _bqp.df_dgradu
    h                       = 1.0./sqrt.(metric.lambda_bar[ParentElems])*ones(1,_bqp.nqp)
    
    #Rotation matrices:
    dq_du                   = dq_du_boundary(nb, tb, nVars)
    du_dq                   = du_dq_boundary(nb, tb, nVars)
    
    #----------------------------------------------------------------
    #Hyperbolic flux. Impose conditions on "u":
    
    #NOTE: The procedure is:
    #1) Rotate u to obtain q.
    #2) From q, define qBC.
    #3) Rotate qBC to obtain uBC.
    #4) Compute flux(uBC) and dfn_I/du_J = dfn/duBC_K * duBC_K/dqBC_L * dqBC_L/dq_M * dq_M/du_J
    
    #Rotate, impose BC and rotate again:
    q                       = uRotate(nb, tb, u) 
    qBC, dqBC_dq            = qBC_(model, BC, t, x, q, ComputeJ)
    uBC                     = qRotate(nb, tb, qBC)
    
    #Matrix duBC_du:
    duBC_du                 = Matrix{Matrix{Float64}}(undef,nVars,nVars)
    if ComputeJ
        alloc!(duBC_du, size(u[1]))
        duBC_dqBC           = du_dq
        for II=1:nVars, JJ=1:nVars, KK=1:nVars, LL=1:nVars
            @mlv    duBC_du[II,JJ]      += duBC_dqBC[II,KK]*dqBC_dq[KK,LL]*dq_du[LL,JJ]
        end
    end

    #Compute hyperbolic flux. Here we employ (n,t) axis, not (x1, x2):
    flux, dflux_duBC, dflux_dgraduBC    = FluxAllocate(nVars, size(u[1]), ComputeJ)
    uBCdep                              = DepVars(model, t, x, uBC, ComputeJ)
    HyperbolicFlux!(model, uBC, uBCdep, ComputeJ, flux, dflux_duBC)
    fn, dfn_duBC, dfn_dgraduBC          = ProjectFlux(flux, dflux_duBC, dflux_dgraduBC, 
                                            nb, ComputeJ)
    
    #Update boundary flux.
    #Note that
    #   df_I/du_J = df_I/duBC_K * duBC_K/du_J:  
    for II=1:nVars
        @mlv    bflux[II]                   += fn[II]
        if ComputeJ
            for JJ=1:nVars, KK=1:nVars
                @mlv    dbflux_du[II,JJ]    += dfn_duBC[II,KK]*duBC_du[KK,JJ]
            end
        end
    end
    
    #Add penalty terms:
    penalty!(model, h, u, uBC, duBC_du, uBCdep, ComputeJ, bflux, dbflux_du)

    #----------------------------------------------------------------
    #Viscous flux. Impose conditions on "f". 
    
    #NOTE: The procedure is:
    #1) Obtain fn_I from u and grad u.
    #2) Rotate fn_I to obtain g_I.
    #3) From g_I, obtain gBC_I.
    #4) Rotate gBC to obtain fnBC.
    #5) Compute dfnBC_I/du_J = dfnBC_I/dgBC_K * dgBC_K/dg_L * dg_L/dfn_M * dfn_M/du_J
    #   and dfnBC_I/d(du_J/dx_k) = dfnBC_I/dgBC_K * dgBC_K/dg_L * dg_L/dfn_M * dfn_M/d(du_J/dx_k)
    
    #Change pointers and reset fluxes:
    dflux_du        = dflux_duBC
    dflux_dgradu    = dflux_dgraduBC
    zero!(flux)
    if ComputeJ
        zero!(dflux_du)
        zero!(dflux_dgradu)
    end
    
    #Compute viscous flux:
    udep    = DepVars(model, t, x, u, ComputeJ)
    nu      = udep[DepVarIndex(model, "nu")][1]
    beta    = udep[DepVarIndex(model, "beta")][1]
    epsilon = udep[DepVarIndex(model, "epsilon")][1]
    ViscousFlux!(model, nu, beta, u, udep, du, ComputeJ, flux, dflux_du, dflux_dgradu)
    epsilonFlux!(model, epsilon, du, ComputeJ, flux, dflux_dgradu)
    
    #Total normal flux (f_{Ij} n_j) through the boundary:
    fn, dfn_du, dfn_dgradu      = ProjectFlux(flux, dflux_du, dflux_dgradu, nb, ComputeJ) 
    
    #Rotate fluxes:
    g                           = uRotate(nb, tb, fn)
    
    #Force zero mass fluxes, zero viscous forces in tangent direction and zero heat flux.
    #Note that v_i tau_ij n_j = v_n tau_n + v_t tau_t = 0, so viscous forces do not 
    #introduce power at the boundary.
    #That is, all components of gBC are zero except the normal viscous force:
    gBC                         = Vector{Matrix{Float64}}(undef, nVars)
    alloc!(gBC, size(u[1]))
    @mlv    gBC[nSpecies+1]     = g[nSpecies+1]     
    #
    #Compute derivatives:
    dgBC_dg                     = Matrix{Matrix{Float64}}(undef, nVars, nVars)
    if ComputeJ
        alloc!(dgBC_dg, size(u[1]))
        @mlv    dgBC_dg[nSpecies+1,nSpecies+1]  = 1.0
    end
    
    #Rotate again:
    fnBC                        = qRotate(nb, tb, gBC)
    
    #Matrix dfnBC/dfn:
    dfnBC_dfn                           = Matrix{Matrix{Float64}}(undef, nVars, nVars)
    if ComputeJ
        dfnBC_dgBC                      = du_dq
        dg_dfn                          = dq_du
        alloc!(dfnBC_dfn, size(u[1]))
        for II=1:nVars, JJ=1:nVars, KK=1:nVars, LL=1:nVars
            @mlv    dfnBC_dfn[II,JJ]    += dfnBC_dgBC[II,KK] * dgBC_dg[KK,LL] *
                                            dg_dfn[LL,JJ]
        end
    end
    
    #Update boundary flux.
    #Note that
    #   dfBC_I/du_J         = dfBC_I/df_K * df_K/du_J:
    #   dfBC_I/d(du_J/dx_j) = dfBC_I/df_K * df_K/d(du_J/dx_j):
    for II=1:nVars
        @mlv    bflux[II]                       += fnBC[II]
        if ComputeJ
            for JJ=1:nVars, KK=1:nVars
                @mlv    dbflux_du[II,JJ]        += dfnBC_dfn[II,KK]*dfn_du[KK,JJ]
            end
            for JJ=1:nVars, KK=1:nVars, kk=1:2
                @mlv    dbflux_dgradu[II,JJ,kk] += dfnBC_dfn[II,KK]*dfn_dgradu[KK,JJ,kk]
            end
        end
    end
    
    return
    
end

#Subsonic inlet. nSpecies + 3 conditions:
#   rhoY_k          = rhoY_k_BC
#   rhov_i          = rhov_i_BC
#   q_j n_j         = 0
#Internal energy is extrapolated.
function bflux!(model::GasModel, BC::SubsonicInlet1, 
    _bqp::TrBintVars, ComputeJ::Bool)
    
    nSpecies                = model.nSpecies
    nVars                   = nSpecies+3
    
    nb                      = _bqp.nb
    tb                      = _bqp.tb
    t                       = _bqp.t
    x                       = _bqp.x
    u                       = _bqp.u
    du                      = _bqp.gradu
    ParentElems             = _bqp.Binteg2D.bmesh.ParentElems
    metric                  = _bqp.Binteg2D.mesh.metric
    bflux                   = _bqp.f
    dbflux_du               = _bqp.df_du
    dbflux_dgradu           = _bqp.df_dgradu
    h                       = 1.0./sqrt.(metric.lambda_bar[ParentElems])*ones(1,_bqp.nqp)
    
    #----------------------------------------------------------------
    #Hyperbolic flux. Impose conditions on "u":
    
    #Evaluate rhoY_k and rhov_i at the boundary:
    uDir                    = BC.fun(t, x)
    
    #Impose BC:
    uBC                     = Vector{Matrix{Float64}}(undef,nVars)
    #Copy density and momentum:
    for II in 1:nSpecies+2
        uBC[II]             = @mlv uDir[II]
    end
    rhoBC                   = sum(uBC[1:nSpecies])
    #Extrapolate internal energy:
    rhoe                    = DepVars(model, t, x, u, ["rhoe"])[1][1]
    #Impose total energy:
    uBC[nSpecies+3]         = @mlv rhoe + 0.5*(uBC[nSpecies+1]*uBC[nSpecies+1] + 
                                uBC[nSpecies+2]*uBC[nSpecies+2])/rhoBC
    
    #Derivatives:
    duBC_du                 = Matrix{Matrix{Float64}}(undef,nVars,nVars)
    if ComputeJ
        alloc!(duBC_du, size(u[1]))
        vx              = DepVars(model, t, x, u, ["vx"])[1][1]
        vy              = DepVars(model, t, x, u, ["vy"])[1][1]
        drhoe_drho      = @mlv 0.5*(vx*vx+vy*vy)
        drhoe_drhovx    = @mlv -vx
        drhoe_drhovy    = @mlv -vy
        drhoe_drhoE     = @mlv 1.0 + 0.0*u[nSpecies+3]
        for ss=1:model.nSpecies
            @mlv duBC_du[nSpecies+3, ss]        = drhoe_drho
        end
        @mlv duBC_du[nSpecies+3, nSpecies+1]    = drhoe_drhovx
        @mlv duBC_du[nSpecies+3, nSpecies+2]    = drhoe_drhovy
        @mlv duBC_du[nSpecies+3, nSpecies+3]    = drhoe_drhoE
    end
    
    #Compute hyperbolic flux. Here we employ (n,t) axis, not (x1, x2):
    flux, dflux_duBC, dflux_dgraduBC    = FluxAllocate(nVars, size(u[1]), ComputeJ)
    uBCdep                              = DepVars(model, t, x, uBC, ComputeJ)
    HyperbolicFlux!(model, uBC, uBCdep, ComputeJ, flux, dflux_duBC)
    fn, dfn_duBC, dfn_dgraduBC          = ProjectFlux(flux, dflux_duBC, dflux_dgraduBC, 
                                            nb, ComputeJ)
    
    #Update boundary flux.
    #Note that
    #   df_I/du_J = df_I/duBC_K * duBC_K/du_J:  
    for II=1:nVars
        @mlv    bflux[II]                   += fn[II]
        if ComputeJ
            for JJ=1:nVars, KK=1:nVars
                @mlv    dbflux_du[II,JJ]    += dfn_duBC[II,KK]*duBC_du[KK,JJ]
            end
        end
    end
    
    #Add penalty terms:
    penalty!(model, h, u, uBC, duBC_du, uBCdep, ComputeJ, bflux, dbflux_du)
    
    #----------------------------------------------------------------
    #Viscous flux. Impose null heat flux:
    
    #Change pointers and reset fluxes:
    dflux_du        = dflux_duBC
    dflux_dgradu    = dflux_dgraduBC
    zero!(flux)
    if ComputeJ
        zero!(dflux_du)
        zero!(dflux_dgradu)
    end
    
    #Compute viscous flux:
    udep    = DepVars(model, t, x, u, ComputeJ)
    epsilon = udep[DepVarIndex(model, "epsilon")][1]
    nu      = udep[DepVarIndex(model, "nu")][1]
    beta    = udep[DepVarIndex(model, "beta")][1]
    ViscousFlux!(model, nu, beta, u, udep, du, ComputeJ, flux, dflux_du, dflux_dgradu)
#     BrennerMassFlux!(model, epsilon, u, udep, du, ComputeJ, flux, dflux_du, dflux_dgradu)
    MassDiffusionFlux!(model, u, udep, du, ComputeJ, flux, dflux_du, dflux_dgradu)
    epsilonFlux!(model, epsilon, du, ComputeJ, flux, dflux_dgradu, IIv=Vector{Int}(1:model.nVars-1))
    
    #Total normal flux (f_{Ij} n_j) through the boundary:
    fn, dfn_du, dfn_dgradu      = ProjectFlux(flux, dflux_du, dflux_dgradu, nb, ComputeJ) 
    
    #Update boundary flux.
    #Note that
    #   dfBC_I/du_J         = dfBC_I/df_K * df_K/du_J:
    #   dfBC_I/d(du_J/dx_j) = dfBC_I/df_K * df_K/d(du_J/dx_j)
    #Note: enforcing no power (v_i tau_ij + q_j) n_j through the boundary is unstable.
    for II=1:nVars
        @mlv    bflux[II]                       += fn[II]
        if ComputeJ
            for JJ=1:nVars
                @mlv    dbflux_du[II,JJ]        += dfn_du[II,JJ]
            end
            for JJ=1:nVars, kk=1:2
                @mlv    dbflux_dgradu[II,JJ,kk] += dfn_dgradu[II,JJ,kk]
            end
        end
    end
    
    return
    
end

#Supersonic inlet. nSpecies + 3 conditions:
#   rhoY_k          = rhoY_k_BC
#   rhov_i          = rhov_i_BC
#   rhoE            = rhoE_BC
function bflux!(model::GasModel, BC::SupersonicInlet1, 
    _bqp::TrBintVars, ComputeJ::Bool)
    
    nSpecies                = model.nSpecies
    nVars                   = nSpecies+3
    
    nb                      = _bqp.nb
    tb                      = _bqp.tb
    t                       = _bqp.t
    x                       = _bqp.x
    u                       = _bqp.u
    du                      = _bqp.gradu
    ParentElems             = _bqp.Binteg2D.bmesh.ParentElems
    metric                  = _bqp.Binteg2D.mesh.metric
    bflux                   = _bqp.f
    dbflux_du               = _bqp.df_du
    dbflux_dgradu           = _bqp.df_dgradu
    h                       = 1.0./sqrt.(metric.lambda_bar[ParentElems])*ones(1,_bqp.nqp)
    
    #----------------------------------------------------------------
    #Hyperbolic flux. Impose conditions on "u":
    
    #Evaluate rhoY_k and rhov_i at the boundary:
    uDir                    = BC.fun(t, x)
    
    #Impose BC:
    uBC                     = Vector{Matrix{Float64}}(undef,nVars)
    #Copy density, momentum and total energy:
    for II in 1:nSpecies+3
        uBC[II]             = @mlv uDir[II]
    end
    
    #Derivatives:
    duBC_du                 = Matrix{Matrix{Float64}}(undef,nVars,nVars)
    if ComputeJ
        alloc!(duBC_du, size(u[1]))
        #derivatives of imposed uBC w.r.t. u are zero.
    end
    
    #Compute hyperbolic flux:
    flux, dflux_duBC, dflux_dgraduBC    = FluxAllocate(nVars, size(u[1]), ComputeJ)
    uBCdep                              = DepVars(model, t, x, uBC, ComputeJ)
    HyperbolicFlux!(model, uBC, uBCdep, ComputeJ, flux, dflux_duBC)
    fn, dfn_duBC, dfn_dgraduBC          = ProjectFlux(flux, dflux_duBC, dflux_dgraduBC, 
                                            nb, ComputeJ)
    
    #Update boundary flux.
    #Note that
    #   df_I/du_J = df_I/duBC_K * duBC_K/du_J:  
    for II=1:nVars
        @mlv    bflux[II]                   += fn[II]
        if ComputeJ
            for JJ=1:nVars, KK=1:nVars
                @mlv    dbflux_du[II,JJ]    += dfn_duBC[II,KK]*duBC_du[KK,JJ]
            end
        end
    end
    
    #Add penalty terms:
    penalty!(model, h, u, uBC, duBC_du, uBCdep, ComputeJ, bflux, dbflux_du)

    #----------------------------------------------------------------
    #Diffusive fluxes:
    
    #Change pointers and reset fluxes:
    dflux_du        = dflux_duBC
    dflux_dgradu    = dflux_dgraduBC
    zero!(flux)
    if ComputeJ
        zero!(dflux_du)
        zero!(dflux_dgradu)
    end
    
    #Compute viscous flux:
    udep            = DepVars(model, t, x, u, ComputeJ)
    epsilon         = udep[DepVarIndex(model, "epsilon")][1]
    nu              = udep[DepVarIndex(model, "nu")][1]
    beta            = udep[DepVarIndex(model, "beta")][1]
    kappa_rho_cv    = udep[DepVarIndex(model, "kappa_rho_cv")][1]
    ViscousFlux!(model, nu, beta, u, udep, du, ComputeJ, flux, dflux_du, dflux_dgradu)
#     BrennerMassFlux!(model, epsilon, u, udep, du, ComputeJ, flux, dflux_du, dflux_dgradu)
    epsilonFlux!(model, epsilon, du, ComputeJ, flux, dflux_dgradu)
    MassDiffusionFlux!(model, u, udep, du, ComputeJ, flux, dflux_du, dflux_dgradu)
    HeatFlux!(model, kappa_rho_cv, u, udep, du, ComputeJ, flux, dflux_du, dflux_dgradu)
    
    #Total normal flux (f_{Ij} n_j) through the boundary:
    fn, dfn_du, dfn_dgradu      = ProjectFlux(flux, dflux_du, dflux_dgradu, nb, ComputeJ) 
    
    #Update boundary flux.
    #Note that
    #   dfBC_I/du_J         = dfBC_I/df_K * df_K/du_J:
    #   dfBC_I/d(du_J/dx_j) = dfBC_I/df_K * df_K/d(du_J/dx_j)
    for II=1:nVars
        @mlv    bflux[II]                       += fn[II]
        if ComputeJ
            for JJ=1:nVars
                @mlv    dbflux_du[II,JJ]        += dfn_du[II,JJ]
            end
            for JJ=1:nVars, kk=1:2
                @mlv    dbflux_dgradu[II,JJ,kk] += dfn_dgradu[II,JJ,kk]
            end
        end
    end
    
    return
    
end

#Subsonic outlet. nSpecies+3 conditons:
#   p               = p_BC
#   (rhoY_k, rhov_i are extrapolated)
#   fmass_(k,j) n_j = 0
#   tau_ij n_j      = 0
function bflux!(model::GasModel0, BC::SubsonicOutlet1,
    _bqp::TrBintVars, ComputeJ::Bool) where {GasModel0<:GasModel}
    
    nSpecies                = model.nSpecies
    nVars                   = nSpecies+3
    
    nb                      = _bqp.nb
    tb                      = _bqp.tb
    t                       = _bqp.t
    x                       = _bqp.x
    u                       = _bqp.u
    du                      = _bqp.gradu
    ParentElems             = _bqp.Binteg2D.bmesh.ParentElems
    metric                  = _bqp.Binteg2D.mesh.metric
    bflux                   = _bqp.f
    dbflux_du               = _bqp.df_du
    dbflux_dgradu           = _bqp.df_dgradu
    h                       = 1.0./sqrt.(metric.lambda_bar[ParentElems])*ones(1,_bqp.nqp)
    
    #----------------------------------------------------------------
    #Hyperbolic flux. Impose conditions on "u":
    
    #Evaluate pressure and gamma at the boundary:
    pDir                    = BC.fun(t, x)[1]
    gamma                   = DepVars(model, t, x, u, ["gamma"])[1][1]
    
    #Impose BC:
    uBC                     = Vector{Matrix{Float64}}(undef,nVars)
    #Copy density and momentum from current values:
    for II in 1:nSpecies+2
        uBC[II]             = @mlv u[II]
    end
    rho                     = sum(u[1:nSpecies])
    #Impose total energy:
    if GasModel0==GasIdeal
        uBC[nSpecies+3]     = @mlv pDir/(gamma-1.0) + 0.5*(u[nSpecies+1]*u[nSpecies+1] + 
                                u[nSpecies+2]*u[nSpecies+2])/rho
    elseif GasModel0==GasFP
        uBC[nSpecies+3]     = @mlv pDir/(gamma-1.0) + model.hfF*u[1] + model.hfP*u[2] + 
                                0.5*(u[nSpecies+1]*u[nSpecies+1] + 
                                    u[nSpecies+2]*u[nSpecies+2])/rho
    else
        error("Unexpected")
    end
    
    #Derivatives:
    duBC_du                 = Matrix{Matrix{Float64}}(undef,nVars,nVars)
    if ComputeJ
        alloc!(duBC_du, size(u[1]))
        for II=1:nSpecies+2
            @mlv duBC_du[II,II]     = 1.0
        end
        II                          = nSpecies+3
        for JJ=1:nSpecies
            @mlv duBC_du[II,JJ]     = - 0.5*(u[nSpecies+1]*u[nSpecies+1] + 
                                            u[nSpecies+2]*u[nSpecies+2])/(rho*rho)
        end
        @mlv duBC_du[II,nSpecies+1] = u[nSpecies+1]/rho
        @mlv duBC_du[II,nSpecies+2] = u[nSpecies+2]/rho
    end
    
    #Compute hyperbolic flux. Here we employ (n,t) axis, not (x1, x2):
    flux, dflux_duBC, dflux_dgraduBC    = FluxAllocate(nVars, size(u[1]), ComputeJ)
    uBCdep                              = DepVars(model, t, x, uBC, ComputeJ)
    HyperbolicFlux!(model, uBC, uBCdep, ComputeJ, flux, dflux_duBC)
    fn, dfn_duBC, dfn_dgraduBC          = ProjectFlux(flux, dflux_duBC, dflux_dgraduBC, 
                                            nb, ComputeJ)
    
    #Update boundary flux.
    #Note that
    #   df_I/du_J = df_I/duBC_K * duBC_K/du_J:  
    for II=1:nVars
        @mlv    bflux[II]                   += fn[II]
        if ComputeJ
            for JJ=1:nVars, KK=1:nVars
                @mlv    dbflux_du[II,JJ]    += dfn_duBC[II,KK]*duBC_du[KK,JJ]
            end
        end
    end
    
    #Add penalty terms:
    penalty!(model, h, u, uBC, duBC_du, uBCdep, ComputeJ, bflux, dbflux_du)
    
    #----------------------------------------------------------------
    #Impose null mass and viscous fluxes:
    
    #Change pointers and reset fluxes:
    dflux_du        = dflux_duBC
    dflux_dgradu    = dflux_dgraduBC
    zero!(flux)
    if ComputeJ
        zero!(dflux_du)
        zero!(dflux_dgradu)
    end
    
    #Compute heat flux:
    udep            = DepVars(model, t, x, u, ComputeJ)
    epsilon         = udep[DepVarIndex(model, "epsilon")][1]
    kappa_rho_cv    = udep[DepVarIndex(model, "kappa_rho_cv")][1]
    HeatFlux!(model, kappa_rho_cv, u, udep, du, ComputeJ, flux, dflux_du, dflux_dgradu)
    epsilonFlux!(model, epsilon, du, ComputeJ, flux, dflux_dgradu, IIv=[model.nVars])
    
    #Total normal flux (f_{Ij} n_j) through the boundary:
    fn, dfn_du, dfn_dgradu      = ProjectFlux(flux, dflux_du, dflux_dgradu, nb, ComputeJ) 
    
    #Update boundary flux:
    #   dfBC_I/du_J         = dfBC_I/df_K * df_K/du_J:
    #   dfBC_I/d(du_J/dx_j) = dfBC_I/df_K * df_K/d(du_J/dx_j):
    for II=1:nVars
        @mlv    bflux[II]                       += fn[II]
        if ComputeJ
            for JJ=1:nVars
                @mlv    dbflux_du[II,JJ]        += dfn_du[II,JJ]
            end
            for JJ=1:nVars, kk=1:2
                @mlv    dbflux_dgradu[II,JJ,kk] += dfn_dgradu[II,JJ,kk]
            end
        end
    end
    
    return
    
end

#Do nothing: nSpecies+3 conditions:
#   fmass_(k,j) n_j = 0
#   tau_ij n_j      = 0
#   q_j n_j         = 0
function bflux!(model::GasModel, BC::Union{DoNothing1,SupersonicOutlet1}, 
    _bqp::TrBintVars, ComputeJ::Bool)
    
    nSpecies                = model.nSpecies
    nVars                   = nSpecies+3
    
    nb                      = _bqp.nb
    tb                      = _bqp.tb
    t                       = _bqp.t
    x                       = _bqp.x
    u                       = _bqp.u
    du                      = _bqp.gradu
    ParentElems             = _bqp.Binteg2D.bmesh.ParentElems
    metric                  = _bqp.Binteg2D.mesh.metric
    bflux                   = _bqp.f
    dbflux_du               = _bqp.df_du
    dbflux_dgradu           = _bqp.df_dgradu
    h                       = 1.0./sqrt.(metric.lambda_bar[ParentElems])*ones(1,_bqp.nqp)
    
    #----------------------------------------------------------------
    #Hyperbolic flux. We extrapolate all the variables, i.e., uBC=u:
    
    uBC                     = u
    duBC_du                 = Matrix{Matrix{Float64}}(undef,nVars,nVars)
    if ComputeJ
        alloc!(duBC_du, size(u[1]))
        for II=1:nVars
            @mlv duBC_du[II,II]         .= 1.0
        end
    end
    
    #Compute hyperbolic flux:
    flux, dflux_duBC, dflux_dgraduBC    = FluxAllocate(nVars, size(u[1]), ComputeJ)
    uBCdep                              = DepVars(model, t, x, uBC, ComputeJ)
    HyperbolicFlux!(model, uBC, uBCdep, ComputeJ, flux, dflux_duBC)
    fn, dfn_duBC, dfn_dgraduBC          = ProjectFlux(flux, dflux_duBC, dflux_dgraduBC, 
                                            nb, ComputeJ)
    
    #Update boundary flux.
    #Note that
    #   df_I/du_J = df_I/duBC_K * duBC_K/du_J:  
    for II=1:nVars
        @mlv    bflux[II]                   += fn[II]
        if ComputeJ
            for JJ=1:nVars
                @mlv    dbflux_du[II,JJ]    += dfn_duBC[II,JJ]
            end
        end
    end
    
    #Add penalty terms:
    penalty!(model, h, u, uBC, duBC_du, uBCdep, ComputeJ, bflux, dbflux_du)
    
    #----------------------------------------------------------------
    #Viscous flux:
    
    #We impose zero viscous fluxes
    
    return
    
end

#SupersonicOutlet1 boundary. 3+nSpecies conditions:
#   m_n             = 0     zero normal velocity
#   fmass_(k,j) n_j = 0     no mass flux, k=1,...,nSpecies
#   tau_nj n_j      = 0
#   v_t             = 0
#   q_j n_j         = 0     no heat flux
function qBC_(model::GasModel, BC::SupersonicOutlet1, 
    t::Float64, x::Vector{MFloat}, q::Vector{MFloat},
    ComputeJ::Bool) where MFloat<:Matrix{Float64}
    
    nSpecies                = model.nSpecies
    nVars                   = model.nVars
    
    #Compute qBC and derivatives:
    qBC                     = Vector{MFloat}(undef,nVars)
    for II in vcat(1:nSpecies, nSpecies+1, nSpecies+3)
        qBC[II]             = @mlv q[II]
    end
    qBC[nSpecies+2]         = @mlv 0.0*q[nSpecies+2]    #rhovt = 0
    dqBC_dq                 = Matrix{MFloat}(undef,nVars,nVars)
    if ComputeJ
        alloc!(dqBC_dq, size(q[1]))
        for II in vcat(1:nSpecies, nSpecies+1, nSpecies+3)
            @mlv    dqBC_dq[II,II]              = 1.0
        end
        @mlv    dqBC_dq[nSpecies+2, nSpecies+2] = 0.0
    end
    
    return qBC, dqBC_dq
    
end
    
#SupersonicOutlet1: nSpecies+3 conditions:
#   fmass_(k,j) n_j = 0
#   tau_nj n_j      = 0
#   v_t             = 0
#   q_j n_j         = 0
function bflux_old!(model::GasModel, BC::SupersonicOutlet1, 
    _bqp::TrBintVars, ComputeJ::Bool)
    
    nSpecies                = model.nSpecies
    nVars                   = nSpecies+3
    
    nb                      = _bqp.nb
    tb                      = _bqp.tb
    t                       = _bqp.t
    x                       = _bqp.x
    u                       = _bqp.u
    du                      = _bqp.gradu
    ParentElems             = _bqp.Binteg2D.bmesh.ParentElems
    metric                  = _bqp.Binteg2D.mesh.metric
    bflux                   = _bqp.f
    dbflux_du               = _bqp.df_du
    dbflux_dgradu           = _bqp.df_dgradu
    h                       = 1.0./sqrt.(metric.lambda_bar[ParentElems])*ones(1,_bqp.nqp)
    
    #Rotation matrices:
    dq_du                   = dq_du_boundary(nb, tb, nVars)
    du_dq                   = du_dq_boundary(nb, tb, nVars)
    
    #----------------------------------------------------------------
    #Hyperbolic flux. We extrapolate all the variables, i.e., uBC=u:
    
    #Rotate, impose tangential velocity and rotate again:
    q                       = uRotate(nb, tb, u) 
    qBC, dqBC_dq            = qBC_(model, BC, t, x, q, ComputeJ)
    uBC                     = qRotate(nb, tb, qBC)
    
    #Matrix duBC_du:
    duBC_du                 = Matrix{Matrix{Float64}}(undef,nVars,nVars)
    if ComputeJ
        alloc!(duBC_du, size(u[1]))
        duBC_dqBC           = du_dq
        for II=1:nVars, JJ=1:nVars, KK=1:nVars, LL=1:nVars
            @mlv    duBC_du[II,JJ]      += duBC_dqBC[II,KK]*dqBC_dq[KK,LL]*dq_du[LL,JJ]
        end
    end
    
    #Compute hyperbolic flux:
    flux, dflux_duBC, dflux_dgraduBC    = FluxAllocate(nVars, size(u[1]), ComputeJ)
    uBCdep                              = DepVars(model, t, x, uBC, ComputeJ)
    HyperbolicFlux!(model, uBC, uBCdep, ComputeJ, flux, dflux_duBC)
    fn, dfn_duBC, dfn_dgraduBC          = ProjectFlux(flux, dflux_duBC, dflux_dgraduBC, 
                                            nb, ComputeJ)
    
    #Update boundary flux.
    #Note that
    #   df_I/du_J = df_I/duBC_K * duBC_K/du_J:  
    for II=1:nVars
        @mlv    bflux[II]                   += fn[II]
#         @mlv    bflux[II]                   += uBC[II]
        if ComputeJ
            for JJ=1:nVars, KK=1:nVars
                @mlv    dbflux_du[II,JJ]    += dfn_duBC[II,KK]*duBC_du[KK,JJ]
#                 @mlv    dbflux_du[II,JJ]    += (II==KK)*duBC_du[KK,JJ]
#                 @mlv    dbflux_du[II,JJ]    += dfn_duBC[II,KK]*(KK==JJ)
            end
        end
    end
    
    #Add penalty terms:
    penalty!(model, h, u, uBC, duBC_du, uBCdep, ComputeJ, bflux, dbflux_du)
    
    #----------------------------------------------------------------
    #Viscous flux. Impose conditions on "f". 
    
    #NOTE: The procedure is:
    #1) Obtain fn_I from u and grad u.
    #2) Rotate fn_I to obtain g_I.
    #3) From g_I, obtain gBC_I.
    #4) Rotate gBC to obtain fnBC.
    #5) Compute dfnBC_I/du_J = dfnBC_I/dgBC_K * dgBC_K/dg_L * dg_L/dfn_M * dfn_M/du_J
    #   and dfnBC_I/d(du_J/dx_k) = dfnBC_I/dgBC_K * dgBC_K/dg_L * dg_L/dfn_M * dfn_M/d(du_J/dx_k)
    
    #Change pointers and reset fluxes:
    dflux_du        = dflux_duBC
    dflux_dgradu    = dflux_dgraduBC
    zero!(flux)
    if ComputeJ
        zero!(dflux_du)
        zero!(dflux_dgradu)
    end
    
    #Compute viscous flux:
    udep    = DepVars(model, t, x, u, ComputeJ)
    nu      = udep[DepVarIndex(model, "nu")][1]
    beta    = udep[DepVarIndex(model, "beta")][1]
    epsilon = udep[DepVarIndex(model, "epsilon")][1]
    ViscousFlux!(model, nu, beta, u, udep, du, ComputeJ, flux, dflux_du, dflux_dgradu)
    epsilonFlux!(model, epsilon, du, ComputeJ, flux, dflux_dgradu)
    
    #Total normal flux (f_{Ij} n_j) through the boundary:
    fn, dfn_du, dfn_dgradu      = ProjectFlux(flux, dflux_du, dflux_dgradu, nb, ComputeJ) 
    
    #Rotate fluxes:
    g                           = uRotate(nb, tb, fn)
    
    #Force zero mass fluxes, zero viscous forces in tangent direction and zero heat flux.
    #Note that v_i tau_ij n_j = v_n tau_n + v_t tau_t = 0, so viscous forces do not 
    #introduce power at the boundary.
    #That is, all components of gBC are zero except the normal viscous force:
    gBC                         = Vector{Matrix{Float64}}(undef, nVars)
    alloc!(gBC, size(u[1]))
    @mlv    gBC[nSpecies+2]     = g[nSpecies+2]     
    #
    #Compute derivatives:
    dgBC_dg                     = Matrix{Matrix{Float64}}(undef, nVars, nVars)
    if ComputeJ
        alloc!(dgBC_dg, size(u[1]))
        @mlv    dgBC_dg[nSpecies+2,nSpecies+2]  = 1.0
    end
    
    #Rotate again:
    fnBC                        = qRotate(nb, tb, gBC)
    
    #Matrix dfnBC/dfn:
    dfnBC_dfn                           = Matrix{Matrix{Float64}}(undef, nVars, nVars)
    if ComputeJ
        dfnBC_dgBC                      = du_dq
        dg_dfn                          = dq_du
        alloc!(dfnBC_dfn, size(u[1]))
        for II=1:nVars, JJ=1:nVars, KK=1:nVars, LL=1:nVars
            @mlv    dfnBC_dfn[II,JJ]    += dfnBC_dgBC[II,KK] * dgBC_dg[KK,LL] *
                                            dg_dfn[LL,JJ]
        end
    end
    
    #Update boundary flux.
    #Note that
    #   dfBC_I/du_J         = dfBC_I/df_K * df_K/du_J:
    #   dfBC_I/d(du_J/dx_j) = dfBC_I/df_K * df_K/d(du_J/dx_j):
    for II=1:nVars
        @mlv    bflux[II]                       += fnBC[II]
        if ComputeJ
            for JJ=1:nVars, KK=1:nVars
                @mlv    dbflux_du[II,JJ]        += dfnBC_dfn[II,KK]*dfn_du[KK,JJ]
            end
            for JJ=1:nVars, KK=1:nVars, kk=1:2
                @mlv    dbflux_dgradu[II,JJ,kk] += dfnBC_dfn[II,KK]*dfn_dgradu[KK,JJ,kk]
            end
        end
    end
    
    return
    
end
