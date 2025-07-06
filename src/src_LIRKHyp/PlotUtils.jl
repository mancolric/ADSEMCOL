#Plot variable "vble"
function PlotNodes3D(solver::Union{SolverData, SolverDataSave}, 
    vble::String;
    color="b", marker=".", markersize=0.5)
    
    return PlotNodes3D(solver, solver.model, vble, color=color, marker=marker, 
        markersize=markersize)
    
end

function PlotNodes3D(solver::Union{SolverData, SolverDataSave}, 
    model::ConstModels, 
    vble::String;
    color="b", marker=".", markersize=0.5)
    
    uNodes      = SolutionAtNodes(solver.u, solver.fes)
    NC          = solver.fes.PSpace.NodesCoords
    xm          = [ view(NC,:,1), view(NC,:,2) ]
    v_plot      = DepVars(model, solver.t, 
                    Vector{<:AMF64}(xm), 
                    Vector{<:AMF64}(uNodes), 
                    [vble])[1][1][:] #First component of [vble], first index
    plot3D(solver.fes.PSpace.NodesCoords[:,1], 
            solver.fes.PSpace.NodesCoords[:,2], 
            v_plot, 
            marker=marker, markersize=markersize, color=color, linestyle="none")
    
    return v_plot
    
end


function PlotNodes(sfun::Function, solver::SolverData, vble::String;
    color="b", marker=".", markersize=0.5)
    
    uNodes      = SolutionAtNodes(solver.u, solver.fes)
    NC          = solver.fes.PSpace.NodesCoords
    xm          = [ view(NC,:,1), view(NC,:,2) ]
    v_plot      = DepVars(solver.model, solver.t, 
                    Vector{<:AMF64}(xm), 
                    Vector{<:AMF64}(uNodes), 
                    [vble])[1][1][:] #First component of [vble], first index
    plotNodes(sfun, v_plot, solver.fes, marker=marker, markersize=markersize, color=color)
    
    return v_plot
    
end

function PlotNodes(sfun::Function, solver::SolverDataSave, model::ConstModels,
    vble::String; color="b", marker=".", markersize=0.5)
    
    uNodes      = SolutionAtNodes(solver.u, solver.fes)
    NC          = solver.fes.PSpace.NodesCoords
    xm          = [ view(NC,:,1), view(NC,:,2) ]
    v_plot      = DepVars(model, solver.t, 
                    Vector{<:AMF64}(xm), 
                    Vector{<:AMF64}(uNodes), 
                    [vble])[1][1][:] #First component of [vble], first index
    plotNodes(sfun, v_plot, solver.fes, marker=marker, markersize=markersize, color=color)
    
    return v_plot
    
end

function PlotContour(solver::Union{SolverDataSave, SolverData}, model::ConstModels,
    vble::String; cmap::String="jet", symmetry::String="none", 
    z1::Float64=NaN, z2::Float64=NaN, Nz::Int=100, delta::Float64=1e-10)
    
    uNodes      = SolutionAtNodes(solver.u, solver.fes)
    NC          = solver.fes.PSpace.NodesCoords
    xm          = [ view(NC,:,1), view(NC,:,2) ]
    v_plot      = DepVars(model, solver.t, 
                    Vector{<:AMF64}(xm), 
                    Vector{<:AMF64}(uNodes), 
                    [vble])[1][1][:] #First component of [vble], first index
    if isnan(z1)
        z1      = minimum(v_plot)-delta
    end
    if isnan(z2)
        z2      = maximum(v_plot)+delta
    end
    PlotContour(v_plot, solver.fes.PSpace, cmap=cmap, symmetry=symmetry, 
        zv=linspace(z1,z2,Nz))
    
    return v_plot
    
end

function PlotContour(solver::SolverData, vble::String; 
    cmap::String="jet", symmetry::String="none", 
    z1::Float64=NaN, z2::Float64=NaN, Nz::Int=100, delta::Float64=1e-10)
    
    PlotContour(solver, solver.model, vble, cmap=cmap, symmetry=symmetry, 
        z1=z1, z2=z2, Nz=Nz, delta=delta)
    
    return v_plot
    
end

function SolutionAtNodes(solver::Union{SolverDataSave, SolverData}, model::ConstModels,
    vble::String)
    
    uNodes      = SolutionAtNodes(solver.u, solver.fes)
    NC          = solver.fes.PSpace.NodesCoords
    xm          = [ view(NC,:,1), view(NC,:,2) ]
    v_plot      = DepVars(model, solver.t, 
                    Vector{<:AMF64}(xm), 
                    Vector{<:AMF64}(uNodes), 
                    [vble])[1][1][:] #First component of [vble], first index
                    
    return v_plot
    
end
