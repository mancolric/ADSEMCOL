#Generate Witherden--Vincent quadrature rules with polyquad.
using DelimitedFiles
include("AuxiliaryFunctions/basic.jl")
include("FemFunctions-v3/FemFunctions.jl")

#Requires to have polyquad installed:
#https://github.com/PyFR/Polyquad

function SearchWVRule(q::Int, N0::Int, Nmax::Int)

    current_dir     = @__DIR__

    # q               = 2
    # N0              = 2     #Starting guess for number of points
    # Nmax            = 4     #Maximum number of nodes
    #Loop until we find a rule:
    N               = N0
    while true
        
        println("q=$q, N=$N")

        #Find quadrature:
        OutputFile1     = "$current_dir/FemFunctions-v3/WitherdenVincentRules/$q-$N-compact.txt"
        OutputFile2     = "$current_dir/FemFunctions-v3/WitherdenVincentRules/$q-$N-expanded.txt"
        OutputFile3     = "$current_dir/FemFunctions-v3/WitherdenVincentRules/expanded/tri/$q-$N.txt"
        run(pipeline(`polyquad find -s tri -q $q -n $N -d double -p -V -w 7200`, stdout=OutputFile1))      #We work with Float64. It does not make sense to find more accurate weights

        #Check if at least one rule has been found:
        if countlines(OutputFile1)>0
            
            #Expand to better format:
            run(pipeline(`polyquad expand -s tri -q $q`, stdin=OutputFile1, stdout=OutputFile2))
            
            #Read quadrature rules:
            aux         = readdlm(OutputFile2)
            rule_ptr    = findall(aux[:,1].=="#")
            NRules      = length(rule_ptr)
            xi0m        = [ Float64.(aux[rule_ptr[ii]+1:rule_ptr[ii]+N,1:2]) for ii=1:NRules ]
            w0v         = [ Float64.(aux[rule_ptr[ii]+1:rule_ptr[ii]+N,3]) for ii=1:NRules ]

            #Evaluate truncation error:
            errv        = zeros(NRules)
            for ii=1:NRules

                #Quadrature in the triangle (0,0) -- (1,0) -- (0,1):
                xim     = @. 0.5 + 0.5*xi0m[ii]
                wv      = @. 0.25*w0v[ii]

                #Basis functions of degree q+1:
                Nm      = TrElem_PKD_Nm(xim, q+1)

                #Evaluate int(Nm_(i,j))-c*delta_(i,0)*delta_(j,0). Here c=int(N_(0,0)) = sqrt(2)*0.5 = 1/sqrt(2)
                errv[ii]        = 0.0
                c               = 1.0/sqrt(2)
                for kk=1:size(Nm,2)
                    errv[ii]    += (dot(wv, Nm[:,kk]) - c*(kk==1))^2
                end        
            end

            #Save rule with smallest error:
            aux         = findmin(errv)
            isel        = aux[2]
            writedlm(OutputFile3, hcat(xi0m[isel], w0v[isel]))

            #Remove auxiliary files:
            rm(OutputFile1)
            rm(OutputFile2)

            #Exit loop:
            break

        end

        #Remove auxiliary file 1:
        rm(OutputFile1)

        #Increase N:
        N           += 1
        if N>Nmax
            break
        end

    end

end
