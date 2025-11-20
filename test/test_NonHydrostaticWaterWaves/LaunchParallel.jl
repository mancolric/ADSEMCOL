#Number of workers and threads is set by SLURM, in particular, by the arguments 
#   --ntasks_per_nodes, --cpus-per-task

#See https://www.tquelch.com/posts/distributed-julia/

include("../../src/AuxiliaryFunctions/basic.jl")
using Distributed
using ClusterManagers
using Base.Threads

#-------------------------------------------------------------------------------
#Load cases:

#Read cases:
ProblemNames    = readdlm("$(@__DIR__)/LaunchParallel.txt")
nProblems       = size(ProblemNames,1)

#For each problem, get number of cases to be solved:
nCases          = zeros(Int, nProblems)
for ii=1:nProblems
    ProblemTable    = readdlm("$(@__DIR__)/LaunchParallel_$(ProblemNames[ii]).txt")
    nCases[ii]      = size(ProblemTable, 1)
end
TotalCases      = sum(nCases)

#Vectors with problem names and local cases:
problems        = Vector{String}(undef, TotalCases)
casesv          = Vector{Int}(undef, TotalCases)
kk              = 0
for ii=1:nProblems, jj=1:nCases[ii]
    global kk       += 1
    problems[kk]    = ProblemNames[ii]
    casesv[kk]      = jj
end

#-------------------------------------------------------------------------------
#Run cases in parallel:

#Number of workers is set from the command line.
# addprocs(2)
# addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))
#Workers are created with the same nb of threads than master.

#Set global variable "computer" in all processors:
for ii in workers()
    fetch(@spawnat ii global ProblemNames=ProblemNames)
#     display("$(ii) $(nthreads()) $(BLAS.get_num_threads())")
end

#Load libraries in all processes:
@everywhere begin
    using Pkg
    Pkg.activate(".")
    for ProblemName in ProblemNames
        include("$(@__DIR__)/$(ProblemName).jl")
    end
end

#Function to launch i-th case in the list for problem ProblemName:
@everywhere function LaunchCase(ProblemName::String, iCase::Int)

    try
        
        #Load table:
        ProblemTable    = readdlm("$(@__DIR__)/LaunchParallel_$(ProblemName).txt")
        if ProblemName=="SolitonRelaxed"
            hp0         = 0.1
            FesOrder    = Int(ProblemTable[iCase,1])
            tf          = Float64(ProblemTable[iCase,2])
            xend        = Float64(ProblemTable[iCase,3])
            SolitonId   = Int(ProblemTable[iCase,4])
            TolS        = Float64(ProblemTable[iCase,5])
            TolT        = Float64(ProblemTable[iCase,6])
            SC          = Int(ProblemTable[iCase,7])
            SolitonRelaxed(hp0, FesOrder, tf=tf, Deltat0=1e-6, xend=xend, 
                SolitonId=SolitonId, 
                SpaceAdapt=true, TolS=TolS, 
                TimeAdapt=true, TolT=TolT, 
                SC=SC, SaveRes=true, Deltat_SaveRes=Inf)
                
        elseif ProblemName=="SolitonExact"
            hp0         = 0.1
            FesOrder    = Int(ProblemTable[iCase,1])
            tf          = Float64(ProblemTable[iCase,2])
            xend        = Float64(ProblemTable[iCase,3])
            alpha       = Float64(ProblemTable[iCase,4])
            TolS        = Float64(ProblemTable[iCase,5])
            TolT        = Float64(ProblemTable[iCase,6])
            SC          = Int(ProblemTable[iCase,7])
            SolitonExact(hp0, FesOrder, tf=tf, Deltat0=1e-6, xend=xend, 
                alpha=alpha, 
                SpaceAdapt=true, TolS=TolS, 
                TimeAdapt=true, TolT=TolT, 
                SC=SC, SaveRes=true, Deltat_SaveRes=Inf)
                
        else
            error("Undefined problem $(ProblemName)")
        end
        
    catch err
    
        PrintToFile("$(@__DIR__)/../../temp/LaunchParallel_$(ProblemName)_$(iCase).err", 
            string(Error2String(err), "\n\n"))
        
    end
    
    return
            
end

#Call pmap:
t0              = time()
# pmap(LaunchCase, problems, casesv)
# rmprocs(workers())
println("$(TotalCases) cases computed in $(time()-t0) seconds")
