#Load problem:
ProblemName     = ARGS[1]
include("$(ProblemName).jl")

try
    if ProblemName=="SolitonRelaxed"
        hp0         = 0.1
        FesOrder    = parse(Int, ARGS[2])
        tf          = parse(Float64, ARGS[3])
        xend        = parse(Float64, ARGS[4])
        SolitonId   = parse(Int, ARGS[5])
        TolS        = parse(Float64, ARGS[6])
        TolT        = parse(Float64, ARGS[7])
        SC          = parse(Int, ARGS[8])
        SolitonRelaxed(hp0, FesOrder, tf=tf, Deltat0=1e-6, xend=xend, 
            SolitonId=SolitonId, 
            SpaceAdapt=true, TolS=TolS, 
            TimeAdapt=true, TolT=TolT, 
            SC=SC, SaveRes=true, Deltat_SaveRes=Inf)
            
    elseif ProblemName=="SolitonExact"
        hp0         = 0.1
        FesOrder    = parse(Int, ARGS[2])
        tf          = parse(Float64, ARGS[3])
        xend        = parse(Float64, ARGS[4])
        alpha       = parse(Int, ARGS[5])
        TolS        = parse(Float64, ARGS[6])
        TolT        = parse(Float64, ARGS[7])
        SC          = parse(Int, ARGS[8])
        SolitonExact(hp0, FesOrder, tf=tf, Deltat0=1e-6, xend=xend, 
            alpha=alpha, 
            SpaceAdapt=true, TolS=TolS, 
            TimeAdapt=true, TolT=TolT, 
            SC=SC, SaveRes=true, Deltat_SaveRes=Inf)
            
    else
        error("Undefined problem $(ProblemName)")
    end
    
catch err
    
    println(err)
#     println(Error2String(err))
            
end
