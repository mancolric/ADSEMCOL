#Load problem:
ProblemName     = ARGS[1]
include("$(ProblemName).jl")

try
    if ProblemName=="DetonationWave"
        tf              = parse(Float64, ARGS[2])
        delta           = parse(Float64, ARGS[3])
        epsilon         = parse(Float64, ARGS[4])
        TolS            = parse(Float64, ARGS[5])
        TolT            = parse(Float64, ARGS[6])
        Deltat_SaveRes  = parse(Float64, ARGS[7])
        SC              = parse(Int, ARGS[8])
        DetonationWave(tf=tf, delta=delta, epsilon=epsilon, 
            TolS=TolS, TolT=TolT, 
            SaveRes=true, Deltat_SaveRes=Deltat_SaveRes, SC=SC)
            
    elseif ProblemName=="ShockFlame_FXP"
        FesOrder        = parse(Int, ARGS[2])
        tf              = parse(Float64, ARGS[3])
        q               = parse(Float64, ARGS[4])
        beta            = parse(Float64, ARGS[5])
        ML              = parse(Float64, ARGS[6])
        epsilon         = parse(Float64, ARGS[7])
        Lx              = parse(Float64, ARGS[8])
        SC              = parse(Int, ARGS[9])
        ShockFlame_FXP(FesOrder=FesOrder, tf=tf, q=q, 
            RThatI=beta, RThatB=beta, ML=ML, epsilon=epsilon, 
            Lx=Lx, 
            SaveRes=true, Deltat_SaveRes=tf/30, SC=SC)
            
    else
        error("Undefined problem $(ProblemName)")
    end
    
catch err
    
    println(Error2String(err))
            
end
