include("PlotResults.jl")

#Initial condition:
function u0fun_Leveque(x::Vector{<:AAF64}, delta::Float64)
    
    rcyl    = @avxt @. sqrt((x[1]-0.5)^2 + x[2]^2)
    rhump   = @avxt @. sqrt((x[1]+0.5)^2 + x[2]^2)
    rcone   = @avxt @. sqrt(x[1]^2 + (x[2]+0.5)^2)
    
    ucyl    = @avxt @. SmoothHeaviside(rcyl-0.25, delta, 1.0, 0.0) * #1 if inside the circle
                        (1.0 - 
                            SmoothHeaviside(abs(x[2])-0.05, delta, 1.0, 0.0)*
                            SmoothHeaviside(x[1]-0.4, delta, 0.0, 1.0) )
                            
    uhump   = @avxt @. cos(1.5*pi*rhump)^3*(rhump<1/3)
    
    ucone   = @avxt @. 0.0*(1.0 - rcone/0.25)*(rcone<0.25)
    
    return [ ucyl + uhump + ucone ]
    
end

#Theoretical solution:
function utheor_Leveque(t::Float64, x::Vector{<:AAF64}, delta::Float64)

    #Rotation matrix:
    omega       = 2*pi
    Q           = [ cos(omega*t)    sin(omega*t)
                    -sin(omega*t)   cos(omega*t) ]
    x0          = [ Q[1,1]*x[1]+Q[1,2]*x[2], 
                    Q[2,1]*x[1]+Q[2,2]*x[2] ]
    return u0fun_Leveque(x0, delta)
end
    
#(Use other function for article):
function Lines_Leveque_old(SC::Int, nb::Int; 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)
    
    FileName    = GetFileName(SC, nb)
    solver      = GetSolver(SC, nb)
    delta       = load(FileName, "delta")
    
    #-----------------------------------------
    #Cylinder:
    
    #Line:
    Nx          = 10000
    x1v         = fill(0.5, Nx)
    x2v         = linspace(-0.4, 0.4, Nx)
    
    #Interpolate solution:
    u_terp,     = SolutionCompute(solver.u, solver.fes, [x1v, x2v])
    utheorv     = utheor_Leveque(0.0, [x1v, x2v], delta)
    
    #Plot:
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
#     plot(x2v, utheorv[1], "k", linewidth=0.5)
    plot(x2v, u_terp[1], "b", linewidth=0.7)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Line1.png", dpi=800, pad_inches=0)
    end
    
    #-----------------------------------------
    #Cylinder + hump:
    
    #Line:
    Nx          = 10000
    x1v         = linspace(-1.0+1e-10, 1.0-1e-10, Nx)
    x2v         = fill(0.0, Nx)
    
    #Interpolate solution:
    u_terp,     = SolutionCompute(solver.u, solver.fes, [x1v, x2v])
    utheorv     = utheor_Leveque(0.0, [x1v, x2v], delta)
    
    #Plot:
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
#     plot(x1v, utheorv[1], "k", linewidth=0.5)
    plot(x1v, u_terp[1], "b", linewidth=0.7)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_Line2.png", dpi=800, pad_inches=0)
    end
    
    uP          = GetPolynomials(solver.u, solver.fes)
    umax, imax  = findmax(uP[1])
    umin, imin  = findmin(uP[1])
    println("umax = ", umax, " at ", solver.fes.PSpace.NodesCoords[imax,:])
    println("umin = ", umin, " at ", solver.fes.PSpace.NodesCoords[imin,:])
    
    return
    
end

function Lines_Leveque(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

#     SCv         = vcat(10095:10096,10082)
#     SCv         = [10120]
#     SCv         = [10145]
    SCv         =[SC]
    
    nSC         = length(SCv)
    colorv      = PyPlotColors("jet2", nSC)
    colorv      = ["b"]
    
    fig1        = PyPlotFigure(w=w+0.7, h=h, bottom=1.0, left=0.8)
    fig2        = PyPlotFigure(w=w, h=h, bottom=1.0, left=0.2)
    fig3        = PyPlotFigure(w=w, h=h, bottom=1.0, left=0.2)
    fig4        = PyPlotFigure(w=w, h=h, bottom=1.0, left=0.2)
    
    FileName    = GetFileName(SCv[1], nb)
    delta       = load(FileName, "delta")
    
    #-----------------------------------------
    #Theoretical and numerical solutions:
    
    #Cylinder:
    Nx          = 1000
    x1v         = fill(0.5, Nx)
    x2v         = linspace(-0.4, 0.4, Nx)
    utheorv     = utheor_Leveque(0.0, [x1v, x2v], delta)
    figure(fig1.number)
    plot(x2v, utheorv[1], "-xr", linewidth=0.5)
    for ii=1:nSC
        SC          = SCv[ii]
        solver      = GetSolver(SC, nb)
        u_terp,     = SolutionCompute(solver.u, solver.fes, [x1v, x2v])
        plot(x2v, u_terp[1], color=colorv[ii], linewidth=0.5)
    end
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(L"x_2")
    ylabel(L"u", rotation=0, labelpad=-0.0)
    if SaveFig
        savefig("$(FigUbi)Leveque_Lines1.png", dpi=800, pad_inches=0)
    end
    
    #Cylinder with zoom:
    Nx          = 100
    x1v         = fill(0.5, Nx)
    x2v         = linspace(0.05-5*delta, 0.05+5*delta, Nx)
    utheorv     = utheor_Leveque(0.0, [x1v, x2v], delta)
    figure(fig2.number)
    plot(x2v, utheorv[1], "-xr", linewidth=0.5)
    for ii=1:nSC
        SC          = SCv[ii]
        solver      = GetSolver(SC, nb)
        u_terp,     = SolutionCompute(solver.u, solver.fes, [x1v, x2v])
        plot(x2v, u_terp[1], color=colorv[ii], linewidth=0.5)
    end
    tick_params(axis="both", which="both", labelsize=TickSize)
    xticks([0.045, 0.05, 0.055], ["0.045", "0.050", "0.055"])
    xlabel(L"x_2")
    legend(["Exact", "Num."], fontsize=TickSize) 
    if SaveFig
        savefig("$(FigUbi)Leveque_Lines2.png", dpi=800, pad_inches=0)
    end
    
    #Cylinder + hump:
    Nx          = 1000
    x1v         = linspace(-1.0+1e-10, 1.0-1e-10, Nx)
    x2v         = fill(0.0, Nx)
    utheorv     = utheor_Leveque(0.0, [x1v, x2v], delta)
    figure(fig3.number)
    plot(x1v, utheorv[1], "-xr", linewidth=0.5)
    for ii=1:nSC
        SC          = SCv[ii]
        solver      = GetSolver(SC, nb)
        u_terp,     = SolutionCompute(solver.u, solver.fes, [x1v, x2v])
        plot(x1v, u_terp[1], color=colorv[ii], linewidth=0.5)
    end
    tick_params(axis="x", which="both", labelsize=TickSize)
#     tick_params(axis="y", left=false, right=false)
    xlabel(L"x_1")
    if SaveFig
        savefig("$(FigUbi)Leveque_Lines3.png", dpi=800, pad_inches=0)
    end
    
    #Cylinder + hump with zoom:
    Nx          = 100
    x1v         = linspace(0.25-5*delta, 0.25+5*delta, Nx)
    x2v         = fill(0.0, Nx)
    utheorv     = utheor_Leveque(0.0, [x1v, x2v], delta)
    figure(fig4.number)
    plot(x1v, utheorv[1], "-xr", linewidth=0.5)
    for ii=1:nSC
        SC          = SCv[ii]
        solver      = GetSolver(SC, nb)
        u_terp,     = SolutionCompute(solver.u, solver.fes, [x1v, x2v])
        plot(x1v, u_terp[1], color=colorv[ii], linewidth=0.5)
    end
    tick_params(axis="both", which="both", labelsize=TickSize)
#     xticks([0.245, 0.25, 0.255], [L"0.245", L"0.250", L"0.255"])
    xticks([0.245, 0.25, 0.255], ["0.245", "0.250", "0.255"])
    xlabel(L"x_1")
    legend(["Exact", "Num."], fontsize=TickSize) 
    if SaveFig
        savefig("$(FigUbi)Leveque_Lines4.png", dpi=800, pad_inches=0)
    end
    
    solver      = GetSolver(SCv[1], nb)
    uP          = GetPolynomials(solver.u, solver.fes)
    umax, imax  = findmax(uP[1])
    umin, imin  = findmin(uP[1])
    println("umax = ", umax, " at ", solver.fes.PSpace.NodesCoords[imax,:])
    println("umin = ", umin, " at ", solver.fes.PSpace.NodesCoords[imin,:])
    
    return
    
end

function Contour_Leveque(SC::Int, nb::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=w, 
    bottom::Float64=0.0, left::Float64=0.0, PlotTime::Bool=false)

    solver      = GetSolver(SC, nb)
    
    if PlotTime
        PyPlotFigure(w=w, h=h+0.7, top=0.7, bottom=0.0, left=0.0, right=0.0)
    else 
        PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
    end
    v_plot      = PlotContour(solver.u[1], solver.fes, cmap="jet")
#     axis("equal")
    tick_params(axis="both", which="both", labelsize=TickSize)
    axis("off")
    if PlotTime
        title(latexstring("t=",sprintf1("%.2f",solver.t)), fontsize=LabelSize)
    end
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(nb).png", dpi=800, pad_inches=0)
    end
    
#     PyPlotFigure(w=w, h=h, bottom=bottom, left=left)
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
    PlotContour(solver.u[1], solver.fes, cmap="jet")
    PlotMesh!(solver.mesh, color="w") 
    tick_params(axis="both", which="both", labelsize=TickSize)
    axis("off")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(nb)_Mesh.png", dpi=800, pad_inches=0)
    end
    
#     PyPlotFigure(w=w, h=h, bottom=bottom, left=left)
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
    PlotContour(solver.u[1], solver.fes, cmap="jet")
    PlotMesh!(solver.mesh, color="w") 
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlim(0.69,0.76)
    ylim(0.04,0.11)
    axis("off")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(nb)_Zoom.png", dpi=800, pad_inches=0)
    end
    
    PyPlotFigure(w=w, h=h, top=0.0, bottom=0.0, left=0.0, right=0.0)
    PlotContour(solver.u[1], solver.fes, cmap="jet")
    PlotMesh!(solver.mesh, color="w") 
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlim(0.40-0.02,0.40+0.05)
    ylim(0.05-0.05,0.05+0.02)
    axis("off")
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_$(nb)_Zoom2.png", dpi=800, pad_inches=0)
    end
    
    println("u_min=", minimum(v_plot), ", umax=", maximum(v_plot), ", err_Linf=", errLq_Leveque(SC,nb,q=Inf))
    
    return
    
end

function Plot3D_Leveque(SC::Int, nb::Int)

    solver      = GetSolver(SC, nb)
    
    figure()
    Plot3D(solver.u[1], solver.fes, cmap="jet")

    return
    
end

function MinMax_Leveque(SC::Int, nb::Int; 
    SaveFig::Bool=false, w::Float64=8.0, h::Float64=6.0)

    FileName    = GetFileName(SC, nb)
    tv          = load(FileName, "tv")
    uminv       = load(FileName, "uminv")
    umaxv       = load(FileName, "umaxv")
    errLqv      = load(FileName, "errLqv")
    errLinfv    = load(FileName, "errLinfv")
    
    PyPlotFigure(w=w, h=h, bottom=1.0, left=1.0, right=0.2, top=0.2)
    semilogy(tv, abs.(umaxv.-1.0), "b")
    semilogy(tv, abs.(uminv), "g")
    tick_params(axis="both", which="both", labelsize=TickSize)
    xlabel(L"t^n", fontsize=LabelSize)
    leg         = vcat(latexstring("|\\max ~ u^{h\\tau,n-}-1|"), 
                        latexstring("|\\min ~ u^{h\\tau,n-}|"))
    ylim(1e-6,1e-2)
    legend(leg, fontsize=TickSize) 
    if SaveFig
        savefig("$(FigUbi)Leveque_Oscillations.png", dpi=800, pad_inches=0)
    end
    
#     figure()
#     plot(tv, abs.(uminv), "g")
#     plot(tv, abs.(umaxv.-1.0), "b")
        
    println("max umax=", maximum(umaxv))
    println("min umin=", maximum(uminv))
    println("max e_Linf=", maximum(errLinfv))
    
    figure()
    semilogy(tv, errLqv, "b")
    semilogy(tv, errLinfv, "k")
    
#     figure()
#     plot(tv, errLqv, "b")
#     figure()
#     plot(tv, errLinfv, "k")
#     
    return
    
end

function OscillationIndicator_Leveque(; 
    SaveFig::Bool=false, w::Float64=8.0, h::Float64=6.0)

    SCv         = 10173:10177
    nb          = 10
    
    nSC         = length(SCv)
    colorv      = PyPlotColors("jet2", nSC)
    
    fig1        = PyPlotFigure(w=w, h=h, bottom=1.0, left=2.0, right=0.1, top=0.2)
    fig2        = PyPlotFigure(w=w, h=h, bottom=1.0, left=1.0, right=0.1, top=0.2)
    leg         = String[]
    for ii=1:nSC
    
        SC          = SCv[ii]
        FileName    = GetFileName(SC, nb)
        tv          = load(FileName, "tv")
        uminv       = load(FileName, "uminv")
        umaxv       = load(FileName, "umaxv")
        errLqv      = load(FileName, "errLqv")
        errLinfv    = load(FileName, "errLinfv")
        Deltav      = @. max(abs.(umaxv.-1.0), abs.(uminv.-0.0))
        solver      = GetSolver(SC,nb)
        TolT        = solver.TolT
        
        figure(fig1.number)
        semilogy(tv, Deltav, color=colorv[ii])
    
        figure(fig2.number)
        semilogy(tv, errLinfv, color=colorv[ii])
    
        leg         = vcat(leg, latexstring("\\mathrm{Tol}_T=", sprintf1("%.2E", TolT)))
        
    end
    
    figure(fig1.number)
    xlabel(L"t", fontsize=LabelSize)
    ylabel(latexstring("\\Delta_{ST}"), fontsize=LabelSize, rotation=0, labelpad=15.0)
#     ylim(1e-6,1e-2)
    legend(leg, fontsize=TickSize, labelspacing=0.1) 
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)Leveque_OscIndicator.png", dpi=800, pad_inches=0)
    end
    
    figure(fig2.number)
    xlabel(L"t", fontsize=LabelSize)
    ylabel(latexstring("e_{ST,\\infty}"), fontsize=LabelSize, rotation=0)
#     ylim(1e-6,1e-2)
    legend(leg, fontsize=TickSize, labelspacing=0.1) 
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
#         savefig("$(FigUbi)Leveque_OscIndicator.png", dpi=800, pad_inches=0)
    end
#    
    return
    
end

function CompareSpaceTimeAdaptation_Leveque(StudyCase::String; 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)
                
    SCvv1   = nothing
    nb      = 1
    if StudyCase=="normal"
        SCvv1   = [ 10031:10034, 
                    vcat(10036:10038,10082),  #No AMA at 10035, 10039
                    vcat(10040:10042,10085) ]
    elseif StudyCase=="N3"
        SCvv1   = [ 10043:10044, 
                    10047:10050, 
                    10051:10053 ]
    elseif StudyCase=="ho"
        SCvv1   = [ 10055:10057, 
                    10059:10062, 
                    10063:10063 ]
    elseif StudyCase=="Linf"
        SCvv1   = [ 10067:10070, 
                    10071:10074, 
                    10075:10078 ] #only 10071 was computed
        nb      = 0
        @warn "nb"
    elseif StudyCase=="normal2" #with factor AMA_TolS
        SCvv1   = [ 10089:10092, 
                    vcat(10094:10096,10082),  #No AMA at 10093, 10097
                    vcat(10098:10100,10085) ]
    elseif StudyCase=="normal3" #with factor AMA_TolS, AMA_AnisoOrder=2
        SCvv1   = [ 10113:10115, 
                    10117:10120,
                    10121:10124 ]
    end
    
    Deltatvv1, hpvv1, 
        errvv1, etavv1,
        tCPUvv1, pvv1       = GetVbles(SCvv1, ["Deltat", "hp", "errL2L2", 
                                "etaL2L2", "tCPU", "FesOrder"], nb=nb)
    
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    leg                     = Legend(pvv1, "p", format="%d")
    loglog(NaN,NaN, "-sk", linewidth=0.5, markersize=3.5)
    loglog(NaN,NaN, "--sk", linewidth=0.5, markersize=3.5)
    leg                     = vcat(leg, L"e_{ST}", latexstring("\\eta_{ST}"))
    PlotXY(Deltatvv1, errvv1, linestyle="solid", marker="s")
    PlotXY(Deltatvv1, etavv1, linestyle="dashed", marker="s")
    grid("on")
    tick_params(axis="both", which="both", labelsize=TickSize)
    legend(leg, fontsize=TickSize)
    xlabel(L"N_t^{-1}")
    if SaveFig
        savefig("$(FigUbi)Leveque_SpaceTimeAdapt1.png", dpi=800, pad_inches=0)
    end
    
    EOCvv1                  = ExpOrderConv(Deltatvv1, errvv1)
    display(EOCvv1)
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    leg                     = Legend(pvv1, "p", format="%d")
    loglog(NaN,NaN, "-sk", linewidth=0.5, markersize=3.5)
    loglog(NaN,NaN, "--sk", linewidth=0.5, markersize=3.5)
    leg                     = vcat(leg, L"e_{ST}", latexstring("\\eta_{ST}"))
    PlotXY(hpvv1, errvv1, linestyle="solid", marker="s")
    PlotXY(hpvv1, etavv1, linestyle="dashed", marker="s")
    grid("on")
    tick_params(axis="both", which="both", labelsize=TickSize)
    legend(leg, fontsize=TickSize)
    xlabel(L"NDOF^{-1/2}")
    if SaveFig
        savefig("$(FigUbi)Leveque_SpaceTimeAdapt2.png", dpi=800, pad_inches=0)
    end
    
    EOCvv1                  = ExpOrderConv(hpvv1, errvv1)
    display(EOCvv1)
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    PlotXY(tCPUvv1, errvv1, linestyle="solid", marker="s")
    grid("on")
    tick_params(axis="both", which="both", labelsize=TickSize)
    legend(leg, fontsize=TickSize)
    xlabel(latexstring(GetString("tCPU")), fontsize=LabelSize)
    ylabel(latexstring(GetString("errL2L2")), rotation=0, fontsize=LabelSize)
    if SaveFig
        savefig("$(FigUbi)Leveque_SpaceTimeAdapt3.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function Errors_Leveque(; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    @warn "not finished"
    
    SC      = 10082
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    PlotErrors!(SC, 1, valid=false)
    xlim(0.05302, 0.05317)
    
    return
    
end

function NElems_Leveque(; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

#     SC      = 10157
#     nb      = 1
    SC      = 10177
    nb      = 10
    
    solver  = GetSolver(SC, nb)
    validv  = ones(Bool, length(solver.tv))
    if true
        @. validv  = solver.validv .== 1 && solver.tv>0.0
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.0, left=1.2)
    plot(solver.tv[validv], solver.nElemsv[validv], "k", linewidth=0.5)
    xlabel(L"t")
#     ylabel(latexstring("\\#\\mathcal{T}^{h,n}"), rotation=0, labelpad=15.0)
    text(-0.25, 7000, latexstring("\\#\\mathcal{T}^{h}"))
    ylim(4000, 8000)
    yticks([4000, 6000, 8000], ["4000", "6000", "8000"])
#     ylim(1e3, 1e4)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)Leveque_NElems.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function NElems2_Leveque(; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    SCv         = 10173:10177
    nb          = 10
    
    nSC         = length(SCv)
    colorv      = PyPlotColors("jet2", nSC)
    
    PyPlotFigure(w=w, h=h, bottom=1.0, left=2.0, right=0.1, top=0.2)
    leg         = String[]
    for ii=1:nSC
    
        SC      = SCv[ii]
        
        solver  = GetSolver(SC, nb)
        validv  = ones(Bool, length(solver.tv))
        if true
            @. validv  = solver.validv .== 1 && solver.tv>0.0
        end
    
        plot(solver.tv[validv], solver.nElemsv[validv], color=colorv[ii], linewidth=0.5)
        
        leg         = vcat(leg, latexstring("\\mathrm{Tol}_T=", sprintf1("%.2E", solver.TolT)))
        
    end
    
#     legend(leg, fontsize=TickSize, labelspacing=0.1) 
    
    xlabel(L"t")
    ylabel(latexstring("\\#\\mathcal{T}^{h}"), rotation=0, labelpad=15.0)
    ylim(3000, 8000)
#     yticks([4000, 6000, 8000], ["4000", "6000", "8000"])
#     ylim(1e3, 1e4)
    tick_params(axis="both", which="both", labelsize=TickSize)
    if SaveFig
        savefig("$(FigUbi)Leveque_NElems.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function CFL_Leveque(; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

    SC      = 10082
    nb      = 1
    
    solver  = GetSolver(SC, nb)
    validv  = ones(Bool, length(solver.tv))
    if true
        @. validv  = solver.validv .== 1
    end
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    semilogy(solver.tv[validv], solver.CFLv[validv], "k", linewidth=0.5)
    xlabel(L"t")
    ylabel(L"CFL", rotation=0)
    if SaveFig
        savefig("$(FigUbi)SC$(SC)_NElems.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function errLq_Leveque(SC::Int, nb::Int; q::Float64=Inf)

    solver      = GetSolver(SC, nb)
    FileName    = GetFileName(SC, nb)
    delta       = load(FileName, "delta")
    
    Integ2D     = TrInt(solver.mesh, 2*(solver.FesOrder+2)+1)
    errLq,      = LqError(Integ2D, solver.u, solver.fes, 
                    FW11((x) -> utheor_Leveque(solver.t, x, delta)), [1.0], q=q)
                    
#     uP          = GetPolynomials(solver.u, solver.fes)
#     errLq,      = LqError(Integ2D, uP, solver.fes.PSpace, 
#                     FW11((x) -> utheor_Leveque(solver.t, x, delta)), [1.0], q=q)

    return errLq
    
end

function GlobalErrorDiscussion(SC::Int; SaveFig::Bool=false, h=10.0, w=10.0)

#     SC          = 10094
    nb          = 1
    
    solver      = GetSolver(SC, nb)
    FileName    = GetFileName(SC, nb)
    
    #Space errors:
    validv      = @. solver.validv == 1
    tv          = solver.tv[validv]
    etaSv       = solver.etaSv[validv]
    
    #Global errors:
    _, etagSv   = _etaS_global(SC,nb)
    
    #Exact errors:
    errLqv      = load(FileName, "errLqv")
#     @. errLqv   = 1.0
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5, top=1.5)
    semilogy(tv, etaSv, ".-b", linewidth=0.7, markersize=1.0)
    semilogy(tv, etagSv, ".-g", linewidth=0.7, markersize=1.0)
    semilogy(unique(tv), errLqv, ".-k", linewidth=0.7, markersize=1.0)
    leg         = ["Estimated", "Accumulated", "Exact"]
    legend(leg)
    xlabel(L"t^n")
    title(string("Values at ",L"t^n"))
    if SaveFig
        savefig("$(FigUbi)GlobalError1.png", dpi=800, pad_inches=0)
    end
    
    #=
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    semilogy(tv, etaSv, ".-b", linewidth=0.7, markersize=1.0)
    semilogy(tv, etagSv, ".-g", linewidth=0.7, markersize=1.0)
    semilogy(unique(tv), errLqv, ".-k", linewidth=0.7, markersize=1.0)
    xlim(0.0,0.01)
    legend(leg)
    xlabel(L"t")
    if SaveFig
        savefig("$(FigUbi)GlobalError2.png", dpi=800, pad_inches=0)
    end
    =#
    
    EtaSv       = TimeL2NormErrors(tv, etaSv)
    EtagSv      = TimeL2NormErrors(tv, etagSv)
    ErrLqv      = TimeL2NormErrors(unique(tv), errLqv)
    
#     PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5, top=1.5)
    semilogy(tv, EtaSv, ".-b", linewidth=0.7, markersize=1.0)
    semilogy(tv, EtagSv, ".-g", linewidth=0.7, markersize=1.0)
    semilogy(unique(tv), ErrLqv, ".-k", linewidth=0.7, markersize=1.0)
    leg         = ["Estimated", "Accumulated", "Exact"]
    legend(leg)
    xlabel(L"t^n")
    title(string(L"L^2", " norm in ", L"[0,t^n]"))
    if SaveFig
        savefig("$(FigUbi)GlobalError3.png", dpi=800, pad_inches=0)
    end
    
    #=
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    semilogy(tv, EtaSv, ".-b", linewidth=0.7, markersize=1.0)
    semilogy(tv, EtagSv, ".-g", linewidth=0.7, markersize=1.0)
    semilogy(unique(tv), ErrLqv, ".-k", linewidth=0.7, markersize=1.0)
    xlim(0.0,0.01)
    legend(leg)
    xlabel(L"t")
    if SaveFig
        savefig("$(FigUbi)GlobalError4.png", dpi=800, pad_inches=0)
    end
    =#
    
end

function CompareGlobalSpaceTimeAdaptation_Leveque(StudyCase::String; 
    SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)
                
    SCvv1   = nothing
    nb      = 1
    if StudyCase=="normal"
        SCvv1   = [ 10031:10034, 
                    vcat(10036:10038,10082),  #No AMA at 10035, 10039
                    vcat(10040:10042,10085) ]
    elseif StudyCase=="N3"
        SCvv1   = [ 10043:10044, 
                    10047:10050, 
                    10051:10053 ]
    elseif StudyCase=="ho"
        SCvv1   = [ 10055:10057, 
                    10059:10062, 
                    10063:10063 ]
    elseif StudyCase=="Linf"
        SCvv1   = [ 10067:10070, 
                    10071:10074, 
                    10075:10078 ] #only 10071 was computed
        nb      = 0
        @warn "nb"
    elseif StudyCase=="normal2" #with factor AMA_TolS
        SCvv1   = [ 10089:10092, 
                    vcat(10094:10096,10082),  #No AMA at 10093, 10097
                    vcat(10098:10100,10085) ]
    end
    
    Deltatvv1, hpvv1, 
        errvv1, etavv1,
        tCPUvv1, pvv1       = GetVbles(SCvv1, ["Deltat", "hp", "errL2L2", 
                                "etaS_global", "tCPU", "FesOrder"], nb=nb)
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    leg                     = Legend(pvv1, "p", format="%d")
    loglog(NaN,NaN, "-sk", linewidth=0.5, markersize=3.5)
    loglog(NaN,NaN, "--sk", linewidth=0.5, markersize=3.5)
    leg                     = vcat(leg, L"e_{ST}", latexstring("Acum. \\eta_{S}"))
    PlotXY(Deltatvv1, errvv1, linestyle="solid", marker="s")
    PlotXY(Deltatvv1, etavv1, linestyle="dashed", marker="s")
    grid("on")
    legend(leg)
    xlabel(L"N_t^{-1}")
    if SaveFig
        savefig("$(FigUbi)Leveque_GlobalSpaceTimeAdapt1.png", dpi=800, pad_inches=0)
    end
    
    EOCvv1                  = ExpOrderConv(Deltatvv1, errvv1)
    display(EOCvv1)
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    leg                     = Legend(pvv1, "p", format="%d")
    loglog(NaN,NaN, "-sk", linewidth=0.5, markersize=3.5)
    loglog(NaN,NaN, "--sk", linewidth=0.5, markersize=3.5)
    leg                     = vcat(leg, L"e_{ST}", latexstring(Accum. "\\eta_{S}"))
    PlotXY(hpvv1, errvv1, linestyle="solid", marker="s")
    PlotXY(hpvv1, etavv1, linestyle="dashed", marker="s")
    grid("on")
    legend(leg)
    xlabel(L"NDOF^{-1/2}")
    if SaveFig
        savefig("$(FigUbi)Leveque_GlobalSpaceTimeAdapt2.png", dpi=800, pad_inches=0)
    end
    
    EOCvv1                  = ExpOrderConv(hpvv1, errvv1)
    display(EOCvv1)
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5)
    PlotXY(tCPUvv1, errvv1, linestyle="solid", marker="s")
    grid("on")
    xlabel(latexstring(GetString("tCPU")))
    ylabel(latexstring(GetString("errL2L2")), rotation=0)
    if SaveFig
        savefig("$(FigUbi)Leveque_GlobalSpaceTimeAdapt3.png", dpi=800, pad_inches=0)
    end
    
    return
    
end

function GlobalErrorDiscussion2(SC::Int; SaveFig::Bool=false, h=10.0, w=10.0)

#     SC          = 10094
    nb          = 1
    
    solver      = GetSolver(SC, nb)
    FileName    = GetFileName(SC, nb)
    
    #Space errors:
    validv      = @. solver.validv == 1
    tv          = solver.tv[validv]
    etaSv       = solver.etaSv[validv]
    _, etagSv   = _etaS_global(SC,nb)
    errLqv      = load(FileName, "errLqv")
#     @. errLqv   = 1.0

    #L2 norm of errors:
    EtaSv       = TimeL2NormErrors(tv, etaSv)
    EtagSv      = TimeL2NormErrors(tv, etagSv)
    ErrLqv      = TimeL2NormErrors(unique(tv), errLqv)
    
    PyPlotFigure(w=w, h=h, bottom=1.5, left=1.5, top=1.5)
    semilogy(tv, etaSv, "b", linewidth=0.7, markersize=1.0)
#     semilogy(tv, EtaSv, "c", linewidth=0.7, markersize=1.0)
    semilogy(tv, etagSv, "r", linewidth=0.7, markersize=1.0)
#     semilogy(tv, EtagSv, "m", linewidth=0.7, markersize=1.0)
    semilogy(unique(tv), errLqv, "k", linewidth=0.7, markersize=1.0)
#     semilogy(unique(tv), ErrLqv, "gray", linewidth=0.7, markersize=1.0)
#     leg         = [latexstring("\\eta_S^n"), string(latexstring("\\eta_S^n"), ", ", L"L^2(0,t^n)", " norm"),
#                    "Accumulated", string("Accumulated", ", ", L"L^2(0,t^n)", " norm"),
#                    "Exact", string("Exact", ", ", L"L^2(0,t^n)", " norm"), ]
    leg         = [latexstring("\\eta_S^n"), 
                   "Accumulated", 
                   "Exact" ]
    legend(leg)
    xlabel(L"t^n")
    if SaveFig
        savefig("$(FigUbi)GlobalError3.png", dpi=800, pad_inches=0)
    end
    
end

function TableResults_Leveque(StudyCase::String)

    #Use StudyCase="normal3"
    
    vbles       = [ "FesOrder",         "TolS_max",        
                    "Deltat",           "hp",  
                    "errLinf_Leveque",  "errL2L2",  "etaST",    
                    "iST",              #"CFLmax",   
                    "tCPU"]
    formats     = [ "%.d",              "%.2E",             
                    "%.2E",             "%.2E",
                    "%.2E",             "%.2E",     "%.2E", 
                    "%.2E",             #"%.2E",     
                    "%.1f"]
                
    header      = string("\$p\$ & \$\\Tol_S^{\\max}=\\Tol_T\$ & \$\\tau\$ & \$h/p\$", 
                        " & \$e_{ST,\\infty}\$ & \$e_{ST}\$ & \$\\Ecal_{ST}\$ ", 
                        " & \$i_{ST}\$",
#                         " & \$\\CFL^{\\max}\$",
                        " & \$t_{\\CPU} [s]\$ \\\\")
    
    
    #------------------------------------------------------------
    
    SCvv0   = NaN
    SCvv1   = NaN
    nb      = 1
    
    if StudyCase=="normal3" #with factor AMA_TolS, AMA_AnisoOrder=2
        SCvv1   = [ 10113:10115, 
                    10117:10120,
                    10121:10124 ]
                    
    elseif StudyCase=="normal4"     #TolS fixed, TolT variable
        SCvv1   = [ 10138:10141, 
                    10142:10145, 
                    10146:10149 ]
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", 
                        repeat("r", length(vbles)), 
                        "} \n", 
                        "\\hline \n", 
                        header, 
                        "\\hline \n")
    for ii=1:length(SCvv1)
        SCv         = SCvv1[ii]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)Leveque_Table.txt", table)
#     write("$(FigUbi)Leveque_Table_$(StudyCase).txt", table)
    
    return
    
end

function TableResults2_Leveque(StudyCase::String)
    
    vbles       = [ "FesOrder",         "TolT",     
                    "hp",               "Deltat",      
                    "errLinf_Leveque",            
                    "errL2L2",          
                    "etaST",            
                    "iST",              "CFLmax",   
                    "errmass_1", 
                    "tCPU"]
    formats     = [ "%.d",              "%.2E",     
                    "%.2E",             "%.2E",             
                    "%.2E",             
                    "%.2E",                
                    "%.2E",            
                    "%.4f",             "%.2E",  
                    "%.2E", 
                    "%.1f"]
                
    header      = string("\$p\$ & \$\\Tol_T\$", 
                        " & \$h/p\$ & \$\\tau\$",
                        " & \$\\Delta_{ST}^{N_t-}\$", 
                        " & \$e_{ST}\$",
                        " & \$\\Ecal_{ST}\$", 
                        " & \$i_{ST}\$",
                        " & \$\\CFL^{\\max}\$",
                        " & \$e^{\\cons}_1\$", 
                        " & \$t_{\\CPU} [s]\$ \\\\")
    
    
    #------------------------------------------------------------
    
    SCvv0   = NaN
    SCvv1   = NaN
    nb      = 1
    
    if StudyCase=="normal_L2" 
        SCvv1   = [ 
#                     10178:10179, 
                    10173:10177, 
#                     10182:10183
                    ] 
        nb      = 10
#         @warn "nb"
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", 
                        repeat("r", length(vbles)), 
                        "} \n", 
                        "\\hline \n", 
                        header, 
                        "\\hline \n")
    for ii=1:length(SCvv1)
        SCv         = SCvv1[ii]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)Leveque_Table2.txt", table)
    
    return
    
end

function TableResults_ProfTimes_Leveque(StudyCase::String)

    vbles       = [ "FesOrder",     "TolS_max", "NDOF",         "tmeanCPU",     "EOC_NDOF_tmeanCPU",
                    "tmeanAMA",     "tmeanJm",  "tmeanSCILU",   "tmeanLS",
                    "LS_mean" ]
                    
    formats     = [ "%d",       "%.2E",      "%.0f",     "%.3f",     "%.1f",
                    "%.3f",     "%.3f",     "%.3f",     "%.3f",
                    "%.1f"  ]
    
    header      = string("\$p\$ & \$\\Tol_S^{\\max}\$ & \$\\NDOF\$",
                        " & \\makecell[cr]{Total}", 
                        " & EOC",
                        " & \\makecell[cr]{Mesh \\\\ adaptation}", 
                        " & \\makecell[cr]{Jacobian \\\\ computation}", 
                        " & \\makecell[cr]{Matrix \\\\ factorization}",
                        " & \\makecell[cr]{Solution of \\\\ linear system}",
                        " & \\makecell[cr]{Linear solver \\\\ iterations}",
                        "\\\\")
                        
    #------------------------------------------------------------
    
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="normal3" #SizeOrder=FesOrder, AnisoOrder=2, with Eq. loop

                        
        SCvv1   = [ 10113:10115, 
                    10117:10120,
                    10121:10124 ]

        nb          = 1
        
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", 
                        repeat("r", length(vbles)), 
                        "} \n", 
                        "\\hline \n", 
                        header, 
                        "\\hline \n")
    for ii=1:length(SCvv1)
        SCv         = SCvv1[ii]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)Leveque_ProfTimes.txt", table)
    
    return
    
end

function TableResults_ProfTimes2_Leveque(StudyCase::String)

    vbles       = [ "FesOrder",     "NDOF",     "tmeanCPU",     
                    "tmeanAMA",     "tmeanJm",  "tmeanSCILU",   "tmeanLS",
                    "LS_mean" ]
                    
    formats     = [ "%d",       "%.0f",     "%.3f",
                    "%.3f",     "%.3f",     "%.3f",     "%.3f",
                    "%.1f"  ]
    
    header      = string("\$p\$ & \$\\NDOF\$ & Total",
                        " & \\makecell[cr]{Mesh \\\\ adaptation}", 
                        " & \\makecell[cr]{Jacobian \\\\ computation}", 
                        " & \\makecell[cr]{Matrix \\\\ factorization}",
                        " & \\makecell[cr]{Solution of \\\\ linear system}",
                        " & \\makecell[cr]{Linear solver \\\\ iterations}",
                        "\\\\")
                        
    #------------------------------------------------------------
    
    SCvv1   = NaN
    nb      = NaN
    
    if StudyCase=="normal3" #SizeOrder=FesOrder, AnisoOrder=2, with Eq. loop
                        
        SCvv1       = [ vcat(10115, 10120, 10124) ]
#         SCvv1       = [ vcat(10114, 10119, 10122) ]

        nb          = 1
        
    elseif StudyCase=="normal_L2" 
        SCvv1   = [ vcat(10181, 10176, 10185) ]
                    #=10178:10179, 
                    10173:10177, 
                    10182:10183
                    ]=# 
        nb      = 10
        
    end
    
    #Save results for IMA and AMA:
    table           = string("\\begin{tabular}{", 
                        repeat("r", length(vbles)), 
                        "} \n", 
                        "\\hline \n", 
                        header, 
                        "\\hline \n")
    for ii=1:length(SCvv1)
        SCv         = SCvv1[ii]
        table_SC    = TableVbles(SCv, vbles, formats, nb=nb)
        table       = string(table, table_SC, "\\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)Leveque_ProfTimes2.txt", table)
    
    return
    
end

function Video_Leveque(SC::Int; SaveFig::Bool=false, w::Float64=8.50, h::Float64=8.50)

#     fig     = PyPlotSubPlots(1, 2, w=w, h=h, top=1.0)
    fig     = PyPlotSubPlots(1, 2, w=w, h=h, left=0.1, right=0.1, bottom=0.1, top=0.7)
    
    kSave   = 0
    for nb=0:1
    
        FileName    = GetFileName(SC, nb)
        solver      = GetSolver(SC, nb)
    
        PyPlot.subplot(1,2,1)
        PyPlot.cla()
        PlotContour(solver.u[1], solver.fes)
#         PlotMesh!(SC, nb, color="w")
        tick_params(axis="both", which="both", labelsize=TickSize)
        axis("off")
#         title(latexstring("t^n=", sprintf1("%.2e", solver.t)),
#             fontsize=10)
        
        PyPlot.subplot(1,2,2)
        PyPlot.cla()
        PlotContour(solver.u[1], solver.fes)
        PlotMesh!(SC, nb, color="w")
        tick_params(axis="both", which="both", labelsize=TickSize)
        axis("off")
        suptitle(latexstring("t^n=", sprintf1("%.2e", solver.t)),
            fontsize=10)
        
        if SaveFig
            savefig("$(FigUbi)SC$(SC)_$(kSave).png", dpi=800, pad_inches=0)
        end
        kSave       += 1
        
        pause(1.0)
        
    end
    
    return
    
end

function TableResults_PrePost_Leveque()

    vbles       = [ "Nt",       "Nt_space_rejected", 
                    "NDOF",     "errL2L2", 
                    "tCPU" ]
    formats     = [ "%d",       "%d", 
                    "%d",       "%.2E", 
                    "%.1f" ]
    
    #------------------------------------------------------------
    
    SCv0    = [ 10186, 10188, 10190 ]
    SCv1    = [ 10187, 10189, 10191 ]
    nb      = 1
    
    #Save results for pre- and post-time-step adaptation approaches:
    table           = string("\\begin{tabular}{", repeat("r", length(vbles)+3), "} \n", 
                        "\\hline \n", 
                        "\$ p\$ ",
                        " & \\makecell[cr]{AMA \\\\ type} ", 
                        " & \$ N_t \$ ",
                        " & \$ N_t^{\\mathrm{rej},S} \$ ",
                        " & NDOF ",
                        " & \$ e_{ST} \$ ",
                        " & \$ t_{\\mathrm{CPU}}[s] \$ ", 
                        " & \\makecell[cr]{\$ t_{\\mathrm{CPU}} \$ \\\\ reduction} ", 
                        "\\\\ \\hline \n")
    
    #Plot table:
    for ii=1:length(SCv0)
        
        #Get just one row:
        table_SC1       = TableVbles([SCv1[ii]], vbles, formats, nb=nb, LastSeparator="\t")
        table_SC0       = TableVbles([SCv0[ii]], vbles, formats, nb=nb, LastSeparator="\t")
        
        #Get CPU time reduction:
        tCPU1           = GetVbles(SCv1[ii], nb, ["tCPU"])[1]
        tCPU0           = GetVbles(SCv0[ii], nb, ["tCPU"])[1]
        p1              = GetVbles(SCv1[ii], nb, ["FesOrder"])[1]
        p0              = GetVbles(SCv0[ii], nb, ["FesOrder"])[1]
        tCPU_red        = (tCPU1-tCPU0)/tCPU1
        table           = string(table, sprintf1("%d", p1), " \t & \t ", 
                                        "Post \t & \t", 
                                        table_SC1, 
                                        "\t & - \t \\\\ \n", 
                                        sprintf1("%d", p0), " \t & \t ", 
                                        "Pre  \t & \t", 
                                        table_SC0, 
                                        "\t & ", sprintf1("%2.1f", tCPU_red*100), "\\%", 
                                        "\t \\\\ \\hline \n")
    end
    table           = string(table, "\\end{tabular}")
    print(table)
    
    write("$(FigUbi)Leveque_PrePost.txt", table)
    
    return
    
end
