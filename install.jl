println("Installing packages")
import Pkg
Pkg.add("Conda")
ENV["PYTHON"]   = ""
Pkg.add("PyCall")
Pkg.build("PyCall")
Pkg.add("PyPlot")
Pkg.build("PyPlot")
# Pkg.add("JLD")
Pkg.add("JLD2")
# Pkg.add("Formatting")
Pkg.add("Format")
Pkg.add("LaTeXStrings")
Pkg.add("IterativeSolvers")
Pkg.add("Krylov")
Pkg.add("Preconditioners")
Pkg.add("LinearMaps")
Pkg.add("Polynomials")
Pkg.add("CuthillMcKee")
Pkg.add("Distributed")
Pkg.add("SharedArrays")
Pkg.add("SparseDiffTools")
Pkg.add("Parameters")
Pkg.add("FastGaussQuadrature")
Pkg.add("IncompleteLU")
Pkg.add("ILUZero")
Pkg.add("AlgebraicMultigrid")
Pkg.add("LoopVectorization")
Pkg.add("Interpolations")
Pkg.add("Optim")
Pkg.add("Humanize")
Pkg.add("Distributions")
Pkg.add("Revise")
Pkg.add("FunctionWrappers")
Pkg.add("JET")
Pkg.add("PProf")
# println("WARN: using Format.jl instead of Formatting.jl")

#Add temp folder:
println("Adding temp folder\n")
CurrentDir  = @__DIR__
mkpath("$(CurrentDir)/test/temp") #no problem if directory already exists

#Create files to save results:
folders     = [ "test_ConvectionDiffusion", "test_CompressibleFlow" ]
for folder in folders

    #Open file and write ResUbi:
    filename    = "$(CurrentDir)/test/$(folder)/ResUbi.jl"
    ResUbiFile  = open(filename, "w")
    DefFolder   = "$(CurrentDir)/results/"
    write(ResUbiFile, "#Path to folder where results are to be saved. Please include final slash (/) character:\n")
    write(ResUbiFile, "ResUbi \t= \"$(DefFolder)\"")
    close(ResUbiFile)
    
    printstyled("File ")
    printstyled(filename, color=:blue, bold=true)
    printstyled(" created. ")
    printstyled("Please change in the file", color=:yellow)
    printstyled(" the directory at which the results for the tests in ")
    printstyled(folder, color=:blue, bold=true)
    printstyled(" are to be saved. Default folder is ")
    printstyled(DefFolder, color=:blue, bold=true)
    printstyled(".\n\n")
    
end

#Create files to save figures and videos:
filename    = "$(CurrentDir)/test/results/ResUbi.jl"
ResUbiFile  = open(filename, "w")
DefFolder   = "$(CurrentDir)/test/figures/"
write(ResUbiFile, "#Path to folder where figures are to be saved. Please include final slash (/) character:\n")
write(ResUbiFile, "FigUbi \t= \"$(DefFolder)\"\n")
write(ResUbiFile, "#Path to folder where videos are to be saved. Please include final slash (/) character:\n")
write(ResUbiFile, "VideosUbi \t= \"$(DefFolder)\"")
close(ResUbiFile)

printstyled("File ")
printstyled(filename, color=:blue, bold=true)
printstyled(" created. ")
printstyled("Please change in the file", color=:yellow)
printstyled(" the directory at which ")
printstyled("figures and videos", color=:blue, bold=true)
printstyled(" are to be saved. Default folder is ")
printstyled(DefFolder, color=:blue, bold=true)
printstyled(".\n\n")
    
printstyled("Installation finished\n", color=:light_green)
