include("basic.jl")

using Distributed

"
    rmallprocs()

Remove all processes, except process 1.
"
function rmallprocs()

    processes   = procs()
    processes   = processes[processes.!=1]
    rmprocs(processes)

end
