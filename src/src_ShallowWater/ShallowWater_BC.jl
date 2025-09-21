function uRotate(model::SWE, nb::Vector{MFloat}, tb::Vector{MFloat}, u::Vector{MFloat}) where MFloat<:AbstractMatrix{Float64}

    nVars       = model.nVars

    q           = Vector{Matrix{Float64}}(undef, nVars)

    #h remains the same:
    q[1]        = @tturbo @. u[1]
    #v_n    = v_x * n_x + v_y * n_y
    q[2]        = @tturbo @. u[2]*nb[1] + u[3]*nb[2]
    #t_n    = v_x * t_x + v_y * t_y
    q[3]        = @tturbo @. u[2]*tb[1] + u[3]*tb[2]
    #b remains the same:
    q[4]        = @tturbo @. u[4]

    return q

end

function qRotate(model::SWE, nb::Vector{MFloat}, tb::Vector{MFloat}, qBC::Vector{MFloat}) where MFloat<:AbstractMatrix{Float64}

    nVars       = model.nVars

    uBC         = Vector{Matrix{Float64}}(undef, nVars)

    #h remains the same:
    uBC[1]      = @tturbo @. qBC[1]
    #v_x    = v_n * n_x + v_t * t_x
    uBC[2]      = @tturbo @. qBC[3]*tb[1] + qBC[2]*nb[1]
    #v_y    = v_n * n_y + v_t * t_y
    uBC[3]      = @tturbo @. qBC[3]*tb[2] + qBC[2]*nb[2]
    #b remains the same:
    uBC[4]      = @tturbo @. qBC[4]
    
    return uBC

end

function _dq_du(model::SWE, nb::Vector{MFloat}, tb::Vector{MFloat}) where MFloat<:AbstractMatrix{Float64}


    #Derivatives of q w.r.t. u:
    dq_du                   = Matrix{Matrix{Float64}}(undef, 4, 4)
    
    dq_du[1,1]              = fill(1.0, size(nb[1]))
    dq_du[1,2]              = fill(0.0, size(nb[1]))
    dq_du[1,3]              = fill(0.0, size(nb[1]))
    dq_du[1,4]              = fill(0.0, size(nb[1]))

    dq_du[2,1]              = fill(0.0, size(nb[1]))
    dq_du[2,2]              = nb[1]
    dq_du[2,3]              = nb[2]
    dq_du[2,4]              = fill(0.0, size(nb[1]))

    dq_du[3,1]              = fill(0.0, size(nb[1]))
    dq_du[3,2]              = tb[1]
    dq_du[3,3]              = tb[2]
    dq_du[3,4]              = fill(0.0, size(nb[1]))
    
    dq_du[4,1]              = fill(0.0, size(nb[1]))
    dq_du[4,2]              = fill(0.0, size(nb[1]))
    dq_du[4,3]              = fill(0.0, size(nb[1]))
    dq_du[4,4]              = fill(0.0, size(nb[1]))

    return dq_du

end


function _duBC_dqBC(model::SWE, nb::Vector{MFloat}, tb::Vector{MFloat}) where MFloat<:AbstractMatrix{Float64}

    #Derivatives of uBC w.r.t. qBC:
    duBC_dqBC           = Matrix{Matrix{Float64}}(undef, 4, 4)

    duBC_dqBC[1,1]      = @tturbo @. 1.0 + nb[1]*0
    duBC_dqBC[1,2]      = 0*nb[1]
    duBC_dqBC[1,3]      = 0*nb[1]
    duBC_dqBC[1,4]      = 0*nb[1]

    duBC_dqBC[2,1]      = 0*nb[1]
    duBC_dqBC[2,2]      = nb[1]
    duBC_dqBC[2,3]      = tb[1]
    duBC_dqBC[2,4]      = 0*nb[1]
    
    duBC_dqBC[3,1]      = 0*nb[1]
    duBC_dqBC[3,2]      = nb[2]
    duBC_dqBC[3,3]      = tb[2]
    duBC_dqBC[3,4]      = 0*nb[1]
    
    duBC_dqBC[4,1]      = 0*nb[1]
    duBC_dqBC[4,2]      = 0*nb[1]
    duBC_dqBC[4,3]      = 0*nb[1]
    duBC_dqBC[4,4]      = 0*nb[1]
    
    return duBC_dqBC

end


#SlipAdiabatic boundary. 3+nSpecies conditions:
function bflux!(model::SWE, BC::SlipAdiabatic,
    _bqp::TrBintVars, ComputeJ::Bool)

    @warn "viscous fluxes are missing"
    
    nVars                   = model.nVars

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
    h                       = 1.0./sqrt.(metric.lambda_bar[ParentElems])*ones(1,_bqp.nqp) #OJO que esto es el tamaño de los elementos

    #----------------------------------------------------------------
    #Hyperbolic flux. Impose conditions on "u":

    #NOTE: The procedure is:
    #1) Rotate u to obtain q.
    #2) From q, define qBC.
    #3) Rotate qBC to obtain uBC.
    #4) Compute flux(uBC) and dfn_I/du_J = dfn/duBC_K * duBC_K/dqBC_L * dqBC_L/dq_M * dq_M/du_J

    #Rotate, impose BC and rotate again:
    q                       = uRotate(model, nb, tb, u)
    qBC                     = Vector{Matrix{Float64}}(undef,4)
    qBC[1]                  = @. 1.0*q[1]
    qBC[2]                  = @. 0.0*q[2]
    qBC[3]                  = @. 1.0*q[3]
    qBC[4]                  = @. 1.0*q[4]
    uBC                     = qRotate(model, nb, tb, qBC)

    #Compute derivatives of uBC wrt u, if required. Note that
    #   duBC/du = duBC/dqBC * dqBC/dq * dq/du:
    duBC_du                 = Matrix{Matrix{Float64}}(undef,nVars,nVars)
    if ComputeJ
    
        #Matrix dqBC/dq is given by present boundary conditions:
        dqBC_dq                     = Matrix{Matrix{Float64}}(undef,nVars,nVars)
        alloc!(dqBC_dq, size(q[1]))
        @tturbo @. dqBC_dq[1,1]     += 1.0
        @tturbo @. dqBC_dq[3,3]     += 1.0
        @tturbo @. dqBC_dq[4,4]     += 1.0
        
        #Matrices duBC/dqBC and dq/du are standard:
        dq_du                   = _dq_du(model,nb,tb)
        duBC_dqBC               = _duBC_dqBC(model,nb,tb)

        #Apply chain rule:
        alloc!(duBC_du, size(u[1]))
        for II=1:nVars, JJ=1:nVars, KK=1:nVars, LL=1:nVars
            @tturbo @. duBC_du[II,JJ]   += duBC_dqBC[II,KK]*dqBC_dq[KK,LL]*dq_du[LL,JJ]
        end
    end

    #Compute hyperbolic flux across the boundary:
    flux, dflux_duBC, dflux_dgraduBC    = FluxAllocate(nVars, size(u[1]), ComputeJ)
    HyperbolicFlux!(model, uBC, ComputeJ, flux, dflux_duBC)
    fn, dfn_duBC, dfn_dgraduBC          = ProjectFlux(flux, dflux_duBC, dflux_dgraduBC,
                                            nb, ComputeJ)

    #Update boundary flux.
    #Note that
    #   df_I/du_J = df_I/duBC_K * duBC_K/du_J:
    for II=1:nVars
        @tturbo @. bflux[II]                    += fn[II]
        if ComputeJ
            for JJ=1:nVars, KK=1:nVars
                @tturbo @. dbflux_du[II,JJ]     += dfn_duBC[II,KK]*duBC_du[KK,JJ]
            end
        end
    end

    return

end


function bflux!(model::SWE, BC::SubsonicInlet1,
                _bqp::TrBintVars, ComputeJ::Bool)

    
    @warn "not reviewed"
    
    nVars                   = model.nVars

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
    h                       = 1.0./sqrt.(metric.lambda_bar[ParentElems])*ones(1,_bqp.nqp) #OJO que esto es el tamaño de los elementos

    #----------------------------------------------------------------
    #Hyperbolic flux. Impose conditions on "u":

    #Evaluate rhoY_k and rhov_i at the boundary:
    uDir                    = BC.fun(t, x)

    #Impose BC:
    uBC                     = Vector{Matrix{Float64}}(undef,nVars)
    #Copy density and momentum:
    #Entrada
    uBC[1]                  = u[1]
    for II in 2:3
        uBC[II]             = @mlv uDir[II-1]
    end


    #Derivatives:
    duBC_du                 = Matrix{Matrix{Float64}}(undef,nVars,nVars)

    if ComputeJ

        # Derivada de uBC respecto de u
        duBC_du                     = Matrix{Matrix{Float64}}(undef, 3,3)

        duBC_du[1,1]                = fill(1.0, size(u[1]))
        duBC_du[1,2]                = fill(0.0, size(u[1]))
        duBC_du[1,3]                = fill(0.0, size(u[1]))

        duBC_du[2,1]                = fill(0.0, size(u[1]))
        duBC_du[2,2]                = fill(0.0, size(u[1]))
        duBC_du[2,3]                = fill(0.0, size(u[1]))

        duBC_du[3,1]                = fill(0.0, size(u[1]))
        duBC_du[3,2]                = fill(0.0, size(u[1]))
        duBC_du[3,3]                = fill(0.0, size(u[1]))

    end

    #Compute hyperbolic flux. Here we employ (n,t) axis, not (x1, x2):
    flux, dflux_duBC, dflux_dgraduBC    = FluxAllocate(nVars, size(u[1]), ComputeJ)
    uBCdep                              = DepVars(model, t, x, uBC, model.DepVars)
    HyperbolicFlux!(model, uBC, ComputeJ, flux, dflux_duBC)
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


end

function bflux!(model::SWE, BC::SubsonicOutlet1,
    _bqp::TrBintVars, ComputeJ::Bool)

    @warn "not reviewed"
    
    nVars                   = model.nVars

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
    h                       = 1.0./sqrt.(metric.lambda_bar[ParentElems])*ones(1,_bqp.nqp) #OJO que esto es el tamaño de los elementos

    #----------------------------------------------------------------
    #Hyperbolic flux. Impose conditions on "u":

    #Evaluate rhoY_k and rhov_i at the boundary:
    uDir                    = BC.fun(t, x)

    #Impose BC:
    uBC                     = Vector{Matrix{Float64}}(undef,nVars)

    #Salida
    uBC[1]                  = uDir[1]
    for II in 2:3
        uBC[II]             = @mlv u[II]
    end

    #Derivatives:
    duBC_du                 = Matrix{Matrix{Float64}}(undef,nVars,nVars)

    if ComputeJ

        # Derivada de uBC respecto de u
        duBC_du                     = Matrix{Matrix{Float64}}(undef, 3,3)

        duBC_du[1,1]                = fill(0.0, size(u[1]))
        duBC_du[1,2]                = fill(0.0, size(u[1]))
        duBC_du[1,3]                = fill(0.0, size(u[1]))

        duBC_du[2,1]                = fill(0.0, size(u[1]))
        duBC_du[2,2]                = fill(1.0, size(u[1]))
        duBC_du[2,3]                = fill(0.0, size(u[1]))

        duBC_du[3,1]                = fill(0.0, size(u[1]))
        duBC_du[3,2]                = fill(0.0, size(u[1]))
        duBC_du[3,3]                = fill(1.0, size(u[1]))

    end

    #Compute hyperbolic flux. Here we employ (n,t) axis, not (x1, x2):
    flux, dflux_duBC, dflux_dgraduBC    = FluxAllocate(nVars, size(u[1]), ComputeJ)
    uBCdep                              = DepVars(model, t, x, uBC, model.DepVars)
    HyperbolicFlux!(model, uBC, ComputeJ, flux, dflux_duBC)
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


end

function bflux!(model::SWE, BC::Union{DoNothing1,SupersonicOutlet1},
                _bqp::TrBintVars, ComputeJ::Bool)

    @warn "not reviewed"
    
    nVars                   = 3

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
#     h                       = 1.0./sqrt.(metric.lambda_bar[ParentElems])*ones(1,_bqp.nqp)

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
    HyperbolicFlux!(model, uBC, ComputeJ, flux, dflux_duBC)
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


    #----------------------------------------------------------------
    #Viscous flux:

    #We impose zero viscous fluxes

    return

end
