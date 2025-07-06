cd(@__DIR__)
include("../test_ConvectionDiffusion/test_ConvectionDiffusion.jl")
include("../test_CompressibleFlow/test_CompressibleFlow.jl")
include("ResUbi.jl")

PyPlot.using3D()

#Default options:
LabelSize   = 10
TickSize    = 8
PyPlot.matplotlib.rc("mathtext",fontset="cm")
PyPlot.matplotlib.rc("font",family="serif")
PyPlot.matplotlib.rc("lines",linewidth=0.5)
PyPlot.matplotlib.rc("font",size=10)

function GetFileName(SC::Int, nb::Int)
    return "$(ResUbi)LIRKHyp_SC$(SC)_$(nb).jld2"
end

function GetSolver(SC::Int, nb::Int)
    FileName    = GetFileName(SC, nb)
    return load(FileName, "solver")
end

function GetVbles(SC::Int, nb::Int, vbles::Vector{String})

    FileName    = GetFileName(SC, nb)
    
    #=
    #Check convergence:
    ConvFlag    = load(FileName, "ConvFlag")
    if ConvFlag<=0
        @warn "SC$(SC) did not converge"
    end
    =#
    
    #Return variables:
    nvbles      = length(vbles)
    xout        = Vector{Any}(undef,nvbles)
    solver      = load(FileName, "solver")
    for ii=1:nvbles
        vble    = vbles[ii]
        if vble=="StudyCase"
            xout[ii]    = load(FileName, "StudyCase")
        elseif vble=="errL2L2"
            xout[ii]    = load(FileName, "errL2L2")
        elseif vble=="etaL2L2"
            xout[ii]    = load(FileName, "etaL2L2")
        elseif vble=="errLqv"
            xout[ii]    = load(FileName, "errLqv")
        elseif vble=="etaST"
            xout[ii]    = solver.etaST
        elseif vble=="errLinf_Leveque"
            if solver.ProblemName!="Leveque"
                error("Current problem is $(solver.ProblemName)")
            end
            umaxv       = load(FileName, "umaxv")
            uminv       = load(FileName, "uminv")
#             xout[ii]    = max( maximum(umaxv)-1.0, 0.0-minimum(uminv) )
            Nt          = length(umaxv)
            xout[ii]    = max(abs(umaxv[Nt]-1.0), abs(uminv[Nt]-0))
        elseif vble=="iST"
            xout[ii]    = solver.etaST/load(FileName, "errL2L2")
        elseif vble=="etaS_global"
            etaS_global,    = _etaS_global(SC,nb)
            xout[ii]        = etaS_global
        elseif vble=="hmean" 
            hmean       = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
            xout[ii]    = hmean
        elseif vble=="hp" #quotient h/p ~ sqrt(1/Ndof)
#             hmean       = 2.0*sqrt(solver.Omega/solver.mesh.nElems/TrElem_Area)
#             xout[ii]    = hmean/solver.FesOrder
#             xout[ii]    = (solver.fes.nDof)^(-1/2)
            xout[ii]    = GetNDof(SC,nb)^(-1/2)
        elseif vble=="tf"
            xout[ii]    = solver.tf
        elseif vble=="tCPU"
            xout[ii]    = solver.tCPU
        elseif vble=="memory"
            xout[ii]    = filesize(FileName)
        elseif vble=="FesOrder"
            xout[ii]    = solver.FesOrder
        elseif vble=="mesh"
            xout[ii]    = solver.mesh
        elseif vble=="nDof"
            xout[ii]    = solver.fes.nDof
        elseif vble=="NDOF"
            xout[ii]    = GetNDof(SC, nb)
        elseif vble=="nElems"
            xout[ii]    = solver.mesh.nElems
        elseif vble=="Deltat"
            xout[ii]    = solver.t/solver.Nt
        elseif vble=="CFLmax"
            xout[ii]    = _CFL_max(SC,nb)
        elseif vble=="CFLmax_corr"
            @warn "deprecated"
            xout[ii]    = _CFL_max_corr(SC,nb)
        elseif vble=="SAType"
            if solver.SpaceAdapt==false || solver.AMA_MaxIter==0
                xout[ii]    = "no"
            elseif solver.AMA_AnisoOrder==0
                xout[ii]    = "IMA"
            else
                xout[ii]    = "AMA"
            end
        elseif vble=="TAType"
            if solver.TimeAdapt
                xout[ii]    = "yes"
            else
                xout[ii]    = "no"
            end
        elseif vble=="TolS_bool"
            if solver.SpaceAdapt && solver.AMA_MaxIter>0
                xout[ii]    = solver.TolS_max
            else
                xout[ii]    = NaN
            end
        elseif vble=="TolT_bool"
            if solver.TimeAdapt
                xout[ii]    = solver.TolT
            else
                xout[ii]    = NaN
            end
        elseif vble=="LS_mean"
            RK          = RK_Coefficients(solver.RKMethod)
#             xout[ii]    = solver.LS_total/(solver.Nt*(RK.stages-1))
            xout[ii]    = solver.LS_total/solver.Nt
        elseif length(vble)>5 && vble[1:5]=="tmean"
            xout[ii]    = getfield(solver, Symbol(string("t", vble[6:length(vble)])))/
                            solver.Nt
        elseif length(vble)>6 && vble[1:6]=="%tmean"
            xout[ii]    = 100*getfield(solver, Symbol(string("t", vble[7:length(vble)])))/
                            solver.tCPU
        elseif vble=="Nt_rejected"
            xout[ii]    = NtRejected(SC,nb)
        elseif vble=="Nt_space_rejected"
            xout[ii]    = NtSpaceRejected(SC,nb)
        elseif length(vble)>8 && vble[1:8]=="errmass_"
            errmassv    = GetMass(SC, nb)-GetMass(SC, 0)
            II          = parse(Int, vble[9])
            xout[ii]    = errmassv[II]
        else
            xout[ii]    = getfield(solver, Symbol(vble))
#             error("Variable $(vble) not supported")
        end
    end
    
    return xout
    
end

function GetString(vble::String)
    
    if vble=="errL2L2"
        return "e_{ST}"
    elseif vble=="etaST"
        return "\\eta_{ST}"
    elseif vble=="hmean" 
        return "h"
    elseif vble=="hp" #quotient h/p ~ sqrt(1/Ndof)
        return "h/p"
    elseif vble=="tf"
        return "t_f"
    elseif vble=="tCPU"
        return "t_{CPU} [s]"
    elseif vble=="memory"
        return "memory"
    elseif vble=="FesOrder"
        return "p"
    elseif vble=="nDof"
        return "N_{DOF}"
    elseif vble=="nElems"
        return "N_{E}"
    elseif vble=="Deltat"
        return "\\Delta t"
    elseif vble=="EOC"
        return "\\mathrm{EOC}"
    elseif vble=="CFLmax"
        return "\\mathrm{CFL}^{\\max{}}"
    else
        error("Variable $(vble) not supported")
    end
    
    return
    
end

function GetVbles(SCvv::Vector{<:AbstractVector{Int}}, vbles::Vector{String}; nb::Int=1)

    nSCvv       = length(SCvv)
    nvbles      = length(vbles)
    
    #Allocate:
    xout            = Vector{Vector{Vector{Float64}}}(undef,nvbles)  #vbles * nSCvv * nSCv
    for ivar=1:nvbles
        xout[ivar]  = Vector{Vector{Float64}}(undef,nSCvv)
    end
    
    #Loop SCvv:
    for ii=1:nSCvv
        
        SCv     = SCvv[ii]
        nSCv    = length(SCv)
        for ivar=1:nvbles
            xout[ivar][ii]  = zeros(nSCv)
        end
        
        #Loop SCv:
        for jj=1:nSCv
            aux     = GetVbles(SCv[jj], nb, vbles)
            for ivar=1:nvbles
                xout[ivar][ii][jj]  = aux[ivar]
            end
        end
        
    end
    
    return xout
    
end

function ExpOrderConv(xvv::Vector{Vector{Float64}}, yvv::Vector{Vector{Float64}})

    nSCvv   = length(xvv)
    EOCvv   = Vector{Vector{Float64}}(undef,nSCvv)
    
    for ii=1:nSCvv
        xv          = xvv[ii]
        yv          = yvv[ii]
        EOCv        = @mlv $diff(log(yv))/$diff(log(xv))
        if xv[1]<xv[length(xv)]
            EOCvv[ii]   = vcat(EOCv, NaN)
        else
            EOCvv[ii]   = vcat(NaN, EOCv)
        end
    end
    
    return EOCvv
    
end

function PlotXY(xvv::Vector{Vector{Float64}}, yvv::Vector{Vector{Float64}};
    linestyle::String="solid", linewidth=0.5, marker::String="", 
    markersize::Float64=3.5, 
    colors::Vector{NTuple{4,Float64}}=PyPlotColors("jet2", length(xvv)),
    scale::String="loglog")

    nSCvv       = length(xvv)
    for ii=1:nSCvv
        
        if scale=="loglog"
            loglog(xvv[ii], yvv[ii], linestyle=linestyle, linewidth=linewidth, marker=marker,
                markersize=markersize, color=colors[ii])
        elseif scale=="semilogx"
            semilogx(xvv[ii], yvv[ii], linestyle=linestyle, linewidth=linewidth, marker=marker,
                markersize=markersize, color=colors[ii])
        else
            error("Unknown plot format $(scale)")
        end
        
    end
    
    return
    
end

function Legend(xvv::Vector{Vector{Float64}}, LegNameLatex::String;
    format::String="%.2e",
    linestyle::String="solid", linewidth=0.5, marker::String="", 
    markersize::Float64=3.5, 
    colors::Vector{NTuple{4,Float64}}=PyPlotColors("jet2", length(xvv)),
    scale::String="loglog")

    nSCvv       = length(xvv)
    leg         = Vector{String}(undef,nSCvv)
    
    for ii=1:nSCvv
        if scale=="loglog"
            loglog([NaN], [NaN], linestyle=linestyle, linewidth=linewidth, marker=marker,
                markersize=markersize, color=colors[ii])
        elseif scale=="semilogx"
            semilogx([NaN], [NaN], linestyle=linestyle, linewidth=linewidth, marker=marker,
                markersize=markersize, color=colors[ii])
        else
            error("Unknown plot format $(scale)")
        end
        leg[ii]     = latexstring(LegNameLatex, "=", sprintf1(format, xvv[ii][1]))
    end
    
    return leg
    
end

function PlotMesh(SC::Int, nb::Int; color::String="k")

    fig     = figure()
    PlotMesh!(SC, nb, color=color)
    axis("equal")
    return fig
    
end

function PlotMesh!(SC::Int, nb::Int; color::String="k")

    mesh,   = GetVbles(SC, nb, ["mesh"])
    PlotMesh!(mesh, color=color)
    return
    
end

function Contour!(SC::Int, nb::Int; cmap::String="jet", PlotMesh::Bool=false, 
    color::String="k", Nz::Int=100, delta::Float64=1e-10)

    solver      = GetSolver(SC, nb)
    
#     figure()
    uP          = GetPolynomials(solver.u, solver.fes)
    umin        = minimum(uP[1])
    umax        = maximum(uP[1])
    PlotContour(solver.u[1], solver.fes, cmap=cmap, 
        zv=linspace(umin-delta, umax+delta, Nz))
    if PlotMesh
        PlotMesh!(solver.mesh, color=color)
    end
    
    return
    
end

function PlotErrors!(SC::Int, nb::Int; valid::Bool=false)

    solver  = GetSolver(SC, nb)
    validv  = ones(Bool, length(solver.tv))
    if valid
        @. validv  = solver.validv .== 1
    end
    semilogy(solver.tv[validv], solver.etaSv[validv], ".-b", linewidth=0.7)
    semilogy(solver.tv[validv], solver.etaTv[validv], ".-g", linewidth=0.7)
    semilogy(solver.tv[validv], solver.etaAv[validv], ".-r", linewidth=0.7)
#     semilogy(solver.tv[validv], solver.nElemsv[validv], ".-k")

    return
    
end

function PlotNElems!(SC::Int, nb::Int; valid::Bool=false)

    solver  = GetSolver(SC, nb)
    validv  = ones(Bool, length(solver.tv))
    if valid
        @. validv  = solver.validv .== 1
    end
    semilogy(solver.tv[validv], solver.nElemsv[validv], ".-k", linewidth=0.7)
    
    return
    
end

function PlotNDof!(SC::Int, nb::Int; valid::Bool=false, color::String=".-b")

    solver  = GetSolver(SC, nb)
    validv  = ones(Bool, length(solver.tv))
    if valid
        @. validv  = solver.validv .== 1
    end
    plot(solver.tv[validv], solver.nDofv[validv], color, linewidth=0.7)
    
    return
    
end

function GetNDof(SC::Int, nb::Int)

    solver      = GetSolver(SC, nb)
    if solver.t==0.0
        return solver.fes.nDof/solver.nVars
    end
    validv      = ones(Bool, length(solver.tv))
    @. validv   = solver.validv .== 1
    nDofv       = solver.nDofv[validv]
    tv          = solver.tv[validv]
    Deltatv     = diff(tv)
    NDof        = sum(nDofv[2:length(nDofv)].*Deltatv)/tv[length(tv)]
    
    return NDof/solver.nVars
    
end

function PlotCFL!(SC::Int, nb::Int; valid::Bool=false, color::String=".-b")

    solver  = GetSolver(SC, nb)
    validv  = ones(Bool, length(solver.tv))
    if valid
        @. validv  = solver.validv .== 1
    end
    semilogy(solver.tv[validv], solver.CFLv[validv], color)
    
end

function PlotDeltat!(SC::Int, nb::Int; valid::Bool=false, color::String=".-b")

    solver  = GetSolver(SC, nb)
    Deltatv = vcat(0, diff(solver.tv))
    validv  = ones(Bool, length(solver.tv))
    if valid
        @. validv  = solver.validv .== 1
    end
    Deltatv[Deltatv.<=0.0]  .= NaN
    semilogy(solver.tv[validv], Deltatv[validv], color)
    
end

function _etaS_global(SC::Int, nb::Int)

    solver              = GetSolver(SC, nb)
    validv              = @. solver.validv == 1
    tv                  = solver.tv[validv]
    etaSv               = solver.etaSv[validv]
    Nt                  = solver.Nt
    Delta_tv            = diff(tv)
    Delta_etaSv         = diff(etaSv)
    
    #Find remeshing and correct Deltat_etaSv:
    aux                 = findall(Delta_tv.==0.0)
    Delta_etaSv[aux]    .= 0.0
    
    #Reconstruct etaSv:
    etaSv_global        = etaSv[1] .+ vcat(0.0, cumsum(Delta_etaSv))

#     figure()
#     semilogy(tv, etaSv_global, ".b")
    
    return sqrt(sum( etaSv_global[2:length(tv)].^2.0 .* diff(tv) )), etaSv_global
    
end

function _CFL_max(SC::Int, nb::Int)

    solver              = GetSolver(SC, nb)
    validv              = @. solver.validv == 1
    CFLv                = solver.CFLv[validv]
    CFLv[isnan.(CFLv)]  .= 0.0
    return maximum(CFLv)
    
end

function _CFL_max_corr(SC::Int, nb::Int)

    solver              = GetSolver(SC, nb)
    validv              = @. solver.validv == 1
    CFLv                = solver.CFLv[validv]
    CFLv[isnan.(CFLv)]  .= 0.0
    
    #The CFL above is based on hElem, with
    #   hElem = 1/Mbar^(1/2)
    #Now we correct it so it is based on 
    #   h1 = 1.0/(Mbar*sf1)^(1/2) = hElem/sf1^(1/2)
    #The new CLF is 
    #   CFL_h1 = a tau/h1 = a tau/hElem hElem/h1 = CFL_hElem sf1^(1/2)
    mesh                = solver.mesh
    sf1max              = maximum(mesh.metric.sf[1])
    #We assume that the CFL max takes place at an element with maximum
    #anisotropy similar to the maximum anisotropy.
    return maximum(CFLv)*sqrt(sf1max)
    
    #The CFL above is based on hElem, with 
    #   hElem   = 1/Mbar^(1/2) = sqrt(K/Khat)
    #Now we correct it so it is based on 
    #   h       = sqrt(K) 
    #The new CFL is
    #   CFL_h   = a tau / h = a tau / hElem * hElem/h = CFL_hElem / sqrt(Khat)
    return maximum(CFLv)/TrElem_Area^(1/2)
    
end

function DisplayTimes(SC::Int, nb::Int)

    solver          = GetSolver(SC, nb)
    DisplayTimes(solver)
    
end

function TimeL2NormErrors(tv::Vector{Float64}, etav::Vector{Float64})

    Nt              = length(tv)
    Etav            = zeros(Nt)
    for ii=2:Nt
#         Etav[ii]    = sqrt( (Etav[ii-1]^2*tv[ii-1] + etav[ii]^2*(tv[ii]-tv[ii-1])) / tv[ii] )
        Etav[ii]    = sqrt( Etav[ii-1]^2 + etav[ii]^2*(tv[ii]-tv[ii-1]) ) 
    end
    Etav[1]         = NaN
    return Etav
    
end

#Remove all intermediate results. Does not remove initial condition nor final result:
function RemoveResults(SC::Int)

    for nb=1:999
        FileName    = GetFileName(SC,nb)
        rm(FileName, force=true)
    end
    
    return 
    
end

function MoveResults(SC::Int; nbv::Array{Int}=Vector{Int}(1:999))

    for nb=nbv
        FileName    = GetFileName(SC,nb)
        ResUbi2     = ""
        if computer=="blasc"
            ResUbi2 = "/data/mancolric/FEM-fluids/Julia-codes/LIRKHyp_heavy/"
        end
        FileName2   = "$(ResUbi2)LIRKHyp_SC$(SC)_$(nb).jld2"
        mv(FileName, FileName2)
    end
    
    return 
    
end

#Write tables:
function TableVbles(SCv::AbstractVector{Int}, vbles::Vector{String}, formats::Vector{String}; 
    nb::Int=1, LastSeparator::String="\t \\\\ \n")

    nSC     = length(SCv)
    nVars   = length(vbles)
    table   = string()
    sep     = vcat(fill("\t & \t",nVars-1), LastSeparator)
    var     = NaN
    for ii=1:nSC
        
        for jj=1:nVars
        
            #Compute numerical value:
            vble    = vbles[jj]
            if length(vble)>3 && vble[1:3]=="EOC"
                if ii==1
                    table   = string(table, "-", sep[jj])
                else
                    T1      = findfirst('_', vble)
                    T2      = findlast('_', vble)
                    vble1   = vble[T1+1:T2-1]
                    vble2   = vble[T2+1:length(vble)]
                    var1    = GetVbles(SCv[ii], nb, [vble1])[1]
                    var1_m1 = GetVbles(SCv[ii-1], nb, [vble1])[1]
                    var2    = GetVbles(SCv[ii], nb, [vble2])[1]
                    var2_m1 = GetVbles(SCv[ii-1], nb, [vble2])[1]
                    var     = log(var2/var2_m1)/log(var1/var1_m1)
                    table   = string(table, sprintf1(formats[jj], var), sep[jj])
                end
            elseif vble[1:1]=="%"
                var     = GetVbles(SCv[ii], nb, [vble])[1]
                table   = string(table, sprintf1(formats[jj], var), "\\%", sep[jj])
            elseif vble=="TolS_bool" || vble=="TolT_bool"
                var     = GetVbles(SCv[ii], nb, [vble])[1]
                if isnan(var)
                    table   = string(table, "No adapt.", sep[jj])
                else
                    table   = string(table, sprintf1(formats[jj], var), sep[jj])
                end
            else
                var     = GetVbles(SCv[ii], nb, [vble])[1]
                table   = string(table, sprintf1(formats[jj], var), sep[jj])
            end
            
        end

    end
    
    return table
    
end

#-------------------------------------------------------------

#Analyze CA:
function PlotCA(SC::Int; SaveFig::Bool=false, w::Float64=8.1, h::Float64=6.0)

#     aux         = load("../temp/LIRKHyp_SC$(SC)_LinearSystem.jld2")
    aux         = load("$(ResUbi)LIRKHyp_SC$(SC)_LinearSystem.jld2")
    A           = aux["A"]
    u0          = aux["u"]
    u           = copy(u0)
    b           = aux["b"]
#     etaST       = aux["etaST"]
    etaST       = 1.0   #so CA is just TolA
    solver      = aux["solver"]
    
    #Find "exact" solution:
    @. u        = u0
    LSOutput    = LS_gmres!(A, u, b, AbsTol=1e-14, Display="final",
                            MaxIter=500)
    u_ex        = copy(u)
    
    #Loop CA:
    N           = 13
#     CAv         = logspace(1e-5, 1e-2, N)
    CAv         = logspace(1e-10, 1e-4, N)
    errv        = zeros(N)
    etav        = zeros(N)
    for ii=1:N
        
        CA          = CAv[ii]
        
        #Find solution for current CA:
        @. u        = u0
        LSOutput    = LS_gmres!(A, u, b, AbsTol=CA*etaST, Display="final", MaxIter=500)
        
        #Find exact error:
        ev          = @. u-u_ex
        ev_views    = GetViews(ev, solver.nVars, solver.fes.nDof)
        Integ2D     = TrInt(solver.mesh, 2*(solver.FesOrder+2)+1)
        errv[ii],   = LqMean(Integ2D, ev_views, solver.fes, solver.nFacts, q=solver.SpaceNorm)
                            
        #Find estimation:
        etav[ii]    = LSOutput[3]

    end
    
    PyPlotFigure(w=w, h=h, bottom=1.2, left=1.0, top=0.1, right=0.1)
    loglog(CAv, etav, "--m")
    loglog(CAv, errv, "m")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(latexstring("\\mathrm{Tol}_A^n"), fontsize=LabelSize)
    legend([latexstring("\\mathcal{E}_A^{h\\tau,n,[i]}"), 
            latexstring("e_A^{h\\tau,n,[i]}")], 
            fontsize=LabelSize, 
            labelspacing=0.1)
    if SaveFig
        savefig("$(FigUbi)TolA_SC$(SC).png", dpi=800, pad_inches=0)
    end
    
    figure()
    plot(CAv, etav./errv, "r")
    
    return
    
end

function PlotCAv(StudyCase::String; SaveFig::Bool=false, w::Float64=8.1, h::Float64=6.0)

    SCvv        = NaN
    
    if StudyCase=="Test1"
    
        SCvv        = [ 20589:20592, 
                        20593:20596,
                        20597:20600 ]
                        
    elseif StudyCase=="SmoothVortex"
    
        SCvv        = [ 50096:50099, 
                        50100:50103,
                        50104:50107 ]
                        
    end
    
    nSCvv       = length(SCvv)
    colorv      = PyPlotColors("jet2", nSCvv)
    leg         = String[]
    
    PyPlotFigure(w=w, h=h, bottom=1.2, left=1.5, top=0.2, right=0.1)
    for ii=1:nSCvv
    
        SCv     = SCvv[ii]
        nSCv    = length(SCv)
        p       = 0
        NDOFv   = zeros(nSCv)
        iAv     = zeros(nSCv)
        
        for jj=1:nSCv
        
            SC          = SCvv[ii][jj]
            
            aux         = load("$(ResUbi)LIRKHyp_SC$(SC)_LinearSystem.jld2")
            A           = aux["A"]
            u0          = aux["u"]
            u           = copy(u0)
            b           = aux["b"]
            etaST       = aux["etaST"]
            solver      = aux["solver"]
            
            #Find "exact" solution:
            @. u        = u0
            LSOutput    = LS_gmres!(A, u, b, AbsTol=1e-14, Display="final",
                                    MaxIter=500)
            u_ex        = copy(u)

            #Find current solution:
            CA          = solver.CA
            @. u        = u0
            LSOutput    = LS_gmres!(A, u, b, AbsTol=CA*etaST, Display="final", MaxIter=500)
    
            #Find exact error:
            ev          = @. u-u_ex
            ev_views    = GetViews(ev, solver.nVars, solver.fes.nDof)
            Integ2D     = TrInt(solver.mesh, 2*(solver.FesOrder+2)+1)
            eA,         = LqMean(Integ2D, ev_views, solver.fes, solver.nFacts, q=solver.SpaceNorm)
                        
            #Find effectivity index:
            etaA        = LSOutput[3]
            iAv[jj]     = etaA/eA
            
            #Find NDOF for current time level:
            Nt          = length(solver.nDofv)
            NDOFv[jj]   = solver.nDofv[Nt]/solver.nVars

            #Find p:
            p           = solver.FesOrder
            
        end
        
        semilogx(NDOFv, iAv, color=colorv[ii], marker="s", markersize=2.5)
        leg     = vcat(leg, latexstring("p=",p))
        
    end
    
    tick_params(axis="both", which="both", labelsize=TickSize)
    if StudyCase=="Test1"
#         xlim(1e3, 2e4)
    elseif StudyCase=="SmoothVortex"
        xlim(1e3, 2e5)
    end
    ylim(1, 18)
    yticks([2.0, 7.0, 12.0, 17.0], ["2", "7", "12", "17"])
    xlabel(latexstring("\\mathrm{NDOF}^n"), fontsize=LabelSize)
    ylabel(latexstring("i_A^{h\\tau,n,[i]}"), fontsize=LabelSize, rotation=0, labelpad=10)
    legend(leg, fontsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)TolA_$(StudyCase).png", dpi=800, pad_inches=0)
    end
    
    return
    
end

#--------------------------------------------------------------
#Mass:

function GetMass(SC::Int, nb::Int)

    #Get solver:
    solver              = GetSolver(SC,nb)
    #this includes:
    #   mesh
    #   fes
    #   nVars
    #   uv
    #   u
    Integ2D                         = TrInt(solver.mesh, 2*(solver.FesOrder+2)+1)
    MII                             = MassMatrix(Integ2D, solver.fes)
    ev                              = zeros(solver.fes.nDof)
    ev[1:solver.fes.PSpace.nDof]    .= 1.0
    
    #Compute (normalized) masses:
    massv           = zeros(solver.nVars)
    for II=1:solver.nVars
        massv[II]   = dot(ev, MII, solver.u[II])/(solver.nFacts[II]*solver.Omega)
    end
    
    return massv
    
end

#--------------------------------------------------------------
#Accepted and rejected time steps:

function NtRejected(SC::Int, nb::Int)

    solver          = GetSolver(SC, nb)
    
    Deltatv         = vcat(0.0, diff(solver.tv)) #Deltat_i = t_i - t_(i-1)
    
#     Nt              = @. $sum( solver.validv==1 && Deltatv!=0.0 )
    Nt_rejected     = @. $sum( solver.validv==0 && Deltatv!=0.0 )
    
    #=
    display(Nt)
    display(solver.Nt)
    display(length(solver.tv))
    display(sum(solver.validv))
    
    figure()
    semilogy(solver.tv, solver.etaSv, "o-b")
#     aux         = @. solver.validv==1
    aux         = @. solver.validv==0 && Deltatv>0.0
    semilogy(solver.tv[aux], solver.etaSv[aux], "+r")
    
    semilogy(solver.tv, solver.etaTv, "o-g")
    aux         = @. solver.validv==0 && Deltatv>0.0
    semilogy(solver.tv[aux], solver.etaTv[aux], "+r")
    
    xlim(0.0, 0.10)
    =#
    
    return Nt_rejected
    
end

function NtSpaceRejected(SC::Int, nb::Int)

    solver          = GetSolver(SC, nb)
    
    Deltatv         = vcat(0.0, diff(solver.tv)) #Deltat_i = t_i - t_(i-1)
    
#     Nt              = @. $sum( solver.validv==1 && Deltatv!=0.0 )
    Nt_rejected     = @. $sum( solver.validv==0 && solver.etaSv>solver.TolS_max && Deltatv!=0.0 )
    
    #=
    display(Nt)
    display(solver.Nt)
    display(length(solver.tv))
    display(sum(solver.validv))
    
    figure()
    semilogy(solver.tv, solver.etaSv, "o-b")
#     aux         = @. solver.validv==1
    aux         = @. solver.validv==0 && Deltatv>0.0
    semilogy(solver.tv[aux], solver.etaSv[aux], "+r")
    
    semilogy(solver.tv, solver.etaTv, "o-g")
    aux         = @. solver.validv==0 && Deltatv>0.0
    semilogy(solver.tv[aux], solver.etaTv[aux], "+r")
    
    xlim(0.0, 0.10)
    =#
    
    return Nt_rejected
    
end
