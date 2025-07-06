include("basic.jl")

mutable struct IMEX_RK

    name        ::String
    AE          ::Matrix{Float64}
    AI          ::Matrix{Float64}
#     bE          ::Vector{Float64}
#     bI          ::Vector{Float64}
    c           ::Vector{Float64}
    bhatE       ::Vector{Float64}
    bhatI       ::Vector{Float64}
    order       ::Int
    stages      ::Int
    GSA         ::Bool
    const_diag  ::Bool  #true if all stages have same diagonal coefficient, except the last one if it is not GSA
    IMEX_RK()   = new()
    
end

function ARK4_3_7L_2_SA1()

    #Implicit method:
    gamma       = 1235/10000
    Am          = zeros(8,8)
    Am[3,2]     = 624185399699/4186980696204
    Am[4,2]     = 1258591069120/10082082980243
    Am[4,3]     = -322722984531/8455138723562
    Am[5,2]     = -436103496990/5971407786587
    Am[5,3]     = -2689175662187/11046760208243
    Am[5,4]     = 4431412449334/12995360898505
    Am[6,2]     = -2207373168298/14430576638973
    Am[6,3]     = 242511121179/3358618340039
    Am[6,4]     = 3145666661981/7780404714551
    Am[6,5]     = 5882073923981/14490790706663
    Am[7,2]     = 0.0
    Am[7,3]     = 9164257142617/17756377923965
    Am[7,4]     = -10812980402763/74029279521829
    Am[7,5]     = 1335994250573/5691609445217
    Am[7,6]     = 2273837961795/8368240463276
    for kk=2:7
        Am[kk,kk]   = gamma
        Am[kk,1]    = Am[kk,2]
    end
    Am[8,:]     .= Am[7,:]

    #Explicit method:
    Bm          = zeros(8,8)
    Bm[2,1]     = 247/1000
    Bm[3,1]     = 247/4000
    Bm[3,2]     = 2694949928731/7487940209513
    Bm[4,1]     = 464650059369/8764239774964
    Bm[4,2]     = 878889893998/2444806327765
    Bm[4,3]     = -952945855348/12294611323341
    Bm[5,1]     = 476636172619/8159180917465
    Bm[5,2]     = -1271469283451/7793814740893
    Bm[5,3]     = -859560642026/4356155882851
    Bm[5,4]     = 1723805262919/4571918432560
    Bm[6,1]     = 6338158500785/11769362343261
    Bm[6,2]     = -4970555480458/10924838743837
    Bm[6,3]     = 3326578051521/2647936831840
    Bm[6,4]     = -880713585975/1841400956686
    Bm[6,5]     = -1428733748635/8843423958496
    Bm[7,1]     = 760814592956/3276306540349
    Bm[7,2]     = 760814592956/3276306540349
    Bm[7,3]     = -47223648122716/6934462133451
    Bm[7,4]     = 71187472546993/9669769126921
    Bm[7,5]     = -13330509492149/9695768672337
    Bm[7,6]     = 11565764226357/8513123442827
    Bm[8,:]     .= Am[8,:]

    return Am, Bm
    
end

function ESDIRK3_2_5L_2_SA()

    Am          = zeros(5,5)
    sq2         = sqrt(2.0)
    for ii=2:5
        Am[ii,ii]   = 9/40
    end
    Am[2,1]     = 9/40
    Am[3,1]     = 9*(1+sq2)/80
    Am[3,2]     = 9*(1+sq2)/80
    Am[4,1]     = (22+15*sq2)/(80*(1+sq2))
    Am[4,2]     = (22+15*sq2)/(80*(1+sq2))
    Am[4,3]     = -7/(40*(1+sq2))
    Am[5,1]     = (2398+1205*sq2)/(2835*(4+3*sq2))
    Am[5,2]     = (2398+1205*sq2)/(2835*(4+3*sq2))
    Am[5,3]     = -2374*(1+2*sq2)/(2835*(5+3*sq2))
    Am[5,4]     = 5827/7560
    
    return Am
    
end

function ESDIRK4_3_7L_2_SA()

    Am          = zeros(7,7)

    for ii=2:7
        Am[ii,ii]   = 1/8
    end
    Am[3,2]     = -39188347878/1513744654945
    Am[4,2]     = 1748874742213/5168247530883
    Am[4,3]     = -1748874742213/5795261096931
    Am[5,2]     = -6429340993097/17896796106705
    Am[5,3]     = 9711656375562/10370074603625
    Am[5,4]     = 1137589605079/3216875020685
    Am[6,2]     = 405169606099/1734380148729
    Am[6,3]     = -264468840649/6105657584947
    Am[6,4]     = 118647369377/6233854714037
    Am[6,5]     = 683008737625/4934655825458
    Am[7,2]     = -5649241495537/14093099002237
    Am[7,3]     = 5718691255176/6089204655961
    Am[7,4]     = 2199600963556/4241893152925
    Am[7,5]     = 8860614275765/11425531467341
    Am[7,6]     = -3696041814078/6641566663007
 
    bv          = zeros(7)
    bv[2]       = 1/4
    bv[3]       = 1200237871921/16391473681546
    bv[4]       = 1/2
    bv[5]       = 395/567
    bv[6]       = 89/126
    bv[7]       = 1.0
    
    for ii=2:7
        Am[ii,1]    = bv[ii]-sum(Am[ii,2:ii])
    end
    
    return Am
    
end

function ESDIRK5_4_8L_2_SA()

    Am          = zeros(8,8)

    for ii=2:8
        Am[ii,ii]   = 1/7
    end
    Am[3,2]     = 1521428834970/8822750406821
    Am[4,2]     = 5338711108027/29869763600956
    Am[4,3]     = 1483184435021/6216373359362
    Am[5,2]     = 2264935805846/12599242299355
    Am[5,3]     = 1330937762090/13140498839569
    Am[5,4]     = -287786842865/17211061626069
    Am[6,2]     = 118352937080/527276862197
    Am[6,3]     = -2960446233093/7419588050389
    Am[6,4]     = -3064256220847/46575910191280
    Am[6,5]     = 6010467311487/7886573591137
    Am[7,2]     = 1134270183919/9703695183946
    Am[7,3]     = 4862384331311/10104465681802
    Am[7,4]     = 1127469817207/2459314315538
    Am[7,5]     = -9518066423555/11243131997224
    Am[7,6]     = -811155580665/7490894181109
    Am[8,2]     = 2162042939093/22873479087181
    Am[8,3]     = -4222515349147/9397994281350
    Am[8,4]     = 3431955516634/4748630552535
    Am[8,5]     = -374165068070/9085231819471
    Am[8,6]     = -1847934966618/8254951855109
    Am[8,7]     = 5186241678079/7861334770480

    bv          = zeros(8)
    bv[2]       = 2/7
    bv[3]       = 5779892736881/11850239716711
    bv[4]       = 150/203
    bv[5]       = 27/46
    bv[6]       = 473/532
    bv[7]       = 30/83
    bv[8]       = 1.0
    
    for ii=2:8
        Am[ii,1]    = bv[ii]-sum(Am[ii,2:ii])
    end
    
    return Am
    
end

function RK_CheckLinearOrder(Am::Array{Float64,2},order::Int)

    #Nb of stages:
    ss      = size(Am,1)
    cv      = sum(Am,dims=2)
    #Loop orders:
    for kk=order
        # println("$kk-th order conditions")
        display(cv.^kk - Am*(kk.*cv.^(kk-1)))
    end

end

"
    RK_coeffs = RK_Coefficients(RKMethod)

Returns IMEX_RK structure ffor the method RKMethod.

Available Runge--Kutta methods:

-IMEXEULER:(1,1,1) IMEX method iin Ascher et al. (1997).

-ASCHER1: (1,2,1) IMEX method iin Ascher et al. (1997).

-ASCHER2: (2,2,2) IMEX method iin Ascher et al. (1997).

-ASCHER3: (4,4,3) IMEX method iin Ascher et al. (1997).

-IMEXSSP3: IMEX method from Conde et al. (2017).

-BPR3: BPR(3,5,3) method from Boscarino et al. (2012).

-KC4A: ARK4(3)6L[2]SA method from Kennedy and Carpenter (2003).

-KC4B: ARK4(3)7L[2]SA1 method from Kennedy and Carpenter (2019).

-CALVO4: LIRK4 method from Calvo et al (2001).

-TSZ3: LIRK-W from Tranquilli et al (2016).

"
function RK_Coefficients(RKMethod::String)

    RKMETHOD        = uppercase(RKMethod)
    RK              = IMEX_RK()
    RK.name         = RKMETHOD
    RK.GSA          = false
    RK.const_diag   = false
    if RKMETHOD=="IMEXEULER"
        RK.AI       = [ 0.0 0.0 ;
                        0.0 1.0 ]
        RK.AE       = [ 0.0 0.0 ;
                        1.0 0.0 ]
        RK.order    = 1
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD == "ASCHER1"
		#From Ascher et al.:
        RK.AI       = [ 0.0 0.0 0.0 ;
                        0.0 1.0 0.0
                        0.0 1.0 0.0 ]
        RK.AE       = [ 0.0 0.0 0.0 ;
                        1.0 0.0 0.0
                        0.0 1.0 0.0 ]
        RK.order    = 1
    elseif RKMETHOD == "ASCHER2"
        #From Ascher's paper:
        aux1        = (2-sqrt(2))/2
        aux2        = 1/2/aux1
        RK.AI       = [ 0.0     0.0         0.0 ;
                        0.0     aux1        0.0 ;
                        0.0     1.0-aux1    aux1 ]
        RK.AE       = [ 0.0     0.0         0.0 ;
                        aux1    0.0         0.0 ;
                        1-aux2  aux2        0.0 ]
        RK.order    = 2
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD == "IMEXSSP3"
		#Gotlieb et al:
        RK.AI       = [ 0.0     0.0     0.0     0.0 ;
                        0.0     1.0     0.0     0.0 ;
                        1/6     -1/3    2/3     0.0 ;
                        1/6     1/6     2/3     0.0 ]
        RK.AE       = [ 0.0     0.0     0.0     0.0 ;
                        1.0     0.0     0.0     0.0 ;
                        0.25    0.25    0.0     0.0 ;
                        1/6     1/6     2/3     0.0 ]
        RK.order    = 3
    elseif RKMETHOD == "EXSSP3"
		#Gotlieb et al:
        RK.AE       = [ 0.0     0.0     0.0     0.0 ;
                        1.0     0.0     0.0     0.0 ;
                        0.25    0.25    0.0     0.0 ;
                        1/6     1/6     2/3     0.0 ]
        RK.AI       = copy(RK.AE)
        RK.order    = 3
	elseif RKMETHOD=="ASCHER3"
		RK.AI		= [ 0.0 	0.0 	0.0 	0.0 	0.0 ;
						0.0 	1/2 	0.0 	0.0 	0.0 ;
						0.0 	1/6 	1/2 	0.0 	0.0 ;
						0.0 	-1/2 	1/2 	1/2 	0.0 ;
						0.0 	3/2 	-3/2 	1/2 	1/2 ]
		RK.AE		= [ 0.0 	0.0 	0.0 	0.0  	0.0;
						1/2 	0.0 	0.0 	0.0  	0.0;
						11/18 	1/18 	0.0 	0.0  	0.0;
						5/6 	-5/6 	1/2 	0.0  	0.0;
                        1/4 	7/4 	3/4 	-7/4  	0.0 ]
#         RK.bhatI    = [ 0.0,    1.0,    0.0,    0.0,    0.0 ] #Does not depend on y^(n+1)
        RK.bhatI    = [ 0.5,    0.0,    0.0,    0.0,    0.5 ] #Not L-stable
#         RK.bhatI    = [ 0.0,    0.0,    3/2,    0.0,    -1/2 ]  #Not L-stable
        RK.bhatE    = copy(RK.bhatI)
        RK.order    = 3
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD=="ASCHER3_E"
		RK.AE		= [ 0.0 	0.0 	0.0 	0.0  	0.0;
						1/2 	0.0 	0.0 	0.0  	0.0;
						11/18 	1/18 	0.0 	0.0  	0.0;
						5/6 	-5/6 	1/2 	0.0  	0.0;
                        1/4 	7/4 	3/4 	-7/4  	0.0 ]
        RK.AI       = copy(RK.AE)
        RK.bhatI    = [ 0.0,    1.0,    0.0,    0.0,    0.0 ]
        RK.bhatE    = copy(RK.bhatI)
        RK.order    = 3
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD=="TEST1_E"
        b1          = 0.7547941406e-1
        b2          = 0.4525669807
        b3          = 0.3591593694
        b4          = 0.9661322979e-1
        b5          = 0.1544737747e-1
        b6          = 0.7336281806e-3
        RK.AE       = [ 0.0 0.0 0.0 0.0 0.0 0.0 0.0
                        1/3 0.0 0.0 0.0 0.0 0.0 0.0
                        0.0 2/3 0.0 0.0 0.0 0.0 0.0
                        0.0 0.0 1.0 0.0 0.0 0.0 0.0
                        0.0 0.0 0.0 0.8 0.0 0.0 0.0
                        0.0 0.0 0.0 0.0 1.0 0.0 0.0 
                        b1  b2  b3  b4  b5  b6  0.0 ]
        RK.AI       = copy(RK.AE)
        RK.bhatI    = RK.AE[7,:]
        RK.bhatE    = copy(RK.bhatI)
        RK.order    = 3
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD=="TEST2_E"
        av  = [ 0.038000000000000
                0.100000000000000
                0.450000000000000
                0.650000000000000
                0.800000000000000
                1.000000000000000
                1.000000000000000
                0.100000000000000
                1.000000000000000
                0.020000000000000
                0.200000000000000 ];

        bv  = [ 0.078156242571258
                0.026668464500122
                0.000370342583619
                0.508019416975715
                0.276789398560433
                0.090247244140683
                0.016033774803200
                0.001959085532335
                0.001571371441140
                0.000076744803796
                0.000102239034082
                0.000005675053007 ]; 
        RK.AE       = [ diagm(-1 => av)     zeros(12,1) 
                        transpose(bv)       0.0         ]
        RK.AI       = copy(RK.AE)
        RK.bhatI    = RK.AE[13,:]
        RK.bhatE    = copy(RK.bhatI)
        RK.order    = 2
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD=="TEST3_E"
        av  = [ 0.082000000000000
                0.130000000000000
                0.750000000000000
                0.510000000000000
                0.590000000000000
                0.400000000000000
                0.340000000000000
                0.300000000000000
                0.130000000000000
                0.060000000000000
                0.050000000000000
                0.040000000000000
                0.100000000000000 ];

        bv  = [ 0.068616748237546
                0.012222604438623
                0.010458669398554
                0.325620288160630
                0.227394916571927
                0.131447938991192
                0.082963215074898
                0.044357759278420
                0.020944614538370
                0.017369644249137
                0.021795547949119
                0.020928594637445
                0.014220981949337
                0.001658476521684 ]; 
        RK.AE       = [ diagm(-1 => av)     zeros(length(bv),1) 
                        transpose(bv)       0.0         ]
        RK.AI       = copy(RK.AE)
        RK.bhatI    = RK.AE[length(bv)+1,:]
        RK.bhatE    = copy(RK.bhatI)
        RK.order    = 2
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD=="TEST_RK14"
        ai_bi = [0                  0.040666901351190
                0.078141050622780   0.173344398566415
                0.104139583599528   0.028278000755788
                0.862070655839388   0.345441344389644
                0.557904001275807   0.201557807105074
                0.657375010196336   0.069843940740099
                0.210879080585046   0.060278020442333
                0.301809357683036   0.028052361631552
                0.200178925860207   0.015734938507497
                0.122870427586218   0.011538983689684
                0.062374580563162   0.012426094378389
                0.062070655839388   0.008987638727357
                0.069357411201343   0.003382512750757
                0.082276055436442   0.000467056963998];  
        M           = size(ai_bi,1)
        av          = ai_bi[2:M,1]
        bv          = ai_bi[:,2]
        RK.AE       = [ diagm(-1 => av)     zeros(M,1) 
                        transpose(bv)       0.0         ]
        RK.AI       = copy(RK.AE)
        RK.bhatI    = RK.AE[M+1,:]
        RK.bhatE    = copy(RK.bhatI)
        RK.order    = 2
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD=="TEST_RK20"
        ai_bi = [0                  0.039492033587098
                0.394117647058824   0.453158516636668
                0.704705882352941   0.295424348439353
                0.963529411764706   0.090067812597009
                0.928235294117647   0.021433097024245
                0.856470588235294   0.004459533258558
                0.750588235294118   0.000889183503304
                0.632941176470588   0.000181096014725
                0.269411764705882   0.000075868688548
                0.070000000000000   0.000108162582941
                0.070000000000000   0.000136911659170
                0.050000000000000   0.000215174168043
                0.020000000000000   0.000744921417995
                0.020000000000000   0.002244068384927
                0.020000000000000   0.005766119216098
                0.020000000000000   0.012255695078178
                0.020000000000000   0.020557161051999
                0.020000000000000   0.025262979590406
                0.020000000000000   0.019982963853801
                0.020000000000000   0.007544353236783];  
        M           = size(ai_bi,1)
        av          = ai_bi[2:M,1]
        bv          = ai_bi[:,2]
        RK.AE       = [ diagm(-1 => av)     zeros(M,1) 
                        transpose(bv)       0.0         ]
        RK.AI       = copy(RK.AE)
        RK.bhatI    = RK.AE[M+1,:]
        RK.bhatE    = copy(RK.bhatI)
        RK.order    = 2
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD=="BPR3"
        #3rd order RK BPR(3,5,3) from Boscarino et al (2012):
        RK.AI       = [ 0.0     0.0     0.0     0.0     0.0
                        0.5     0.5     0.0     0.0     0.0 
                        5/18    -1/9    0.5     0.0     0.0 
                        0.5     0.0     0.0     0.5     0.0 
                        0.25    0.0     0.75    -0.5    0.5 ]
        RK.AE       = [ 0.0     0.0     0.0     0.0     0.0
                        1.0     0.0     0.0     0.0     0.0 
                        4/9     2/9     0.0     0.0     0.0 
                        0.25    0.0     0.75    0.0     0.0 
                        0.25    0.0     0.75    0.0     0.0 ]   
#         RK.bhatI    = [ 0.5     0.0     0.0     0.0     0.5 ] #Not L-stable
#         RK.bhatE    = [ 0.5     0.0     0.0     0.0     0.5 ]
        #Solution of stage 3: second order accurate, L-stable, takes into account info
        #at t^(n+1), good stability region with EINRK:
        RK.bhatI    = [ 5/18    -1/9    0.5     0.0     0.0 ]   
        RK.bhatE    = [ 4/9     2/9     0.0     0.0     0.0 ]
        RK.order    = 3  
        RK.GSA          = true
        RK.const_diag   = true
    # elseif RKMETHOD=="DIRK3"
        # #From wikipedia. No references found:
        # Am 			= [ 0.0 	0.0 	0.0 	0.0 	0.0 ;
        #                 0.0 	1/2 	0.0 	0.0 	0.0 ;
        #                 0.0 	1/6 	1/2 	0.0 	0.0 ;
        #                 0.0 	-1/2 	1/2 	1/2 	0.0 ;
        #                 0.0 	3/2 	-3/2 	1/2 	1/2 ]
        # Bm 			= Am
    elseif RKMETHOD=="KC4A"
		#4-th order Kennedy-Carpenter ARK4(3)6L[2]SA from ANM2013:
        RK.AI       = [ 0                               0                               0                               0                               0           0       0
                        0.25                            0.25                            0                               0                               0           0       0
                        8611/62500                      -1743/31250                     0.25                            0                               0           0       0
                        5012029/34652500                -654441/2922500                 174375/388108                   0.25                            0           0       0
                        15267082809/155376265600        -71443401/120774400             730878875/902184768             2285395/8070912                 0.25        0       0
                        82889/524892                    0                               15625/83664                     69875/102672                    -2260/8211  0.25    0
                        82889/524892                    0                               15625/83664                     69875/102672                    -2260/8211  0.25    0   ]
        RK.AE       = [ 0                               0                               0                               0                               0           0       0
                        0.5                             0                               0                               0                               0           0       0
                        13861/62500                     6889/62500                      0                               0                               0           0       0
                        -116923316275/2393684061468     -2731218467317/15368042101831   9408046702089/11113171139209    0                               0           0       0
                        -451086348788/2902428689909    -2682348792572/7519795681897     12662868775082/11960479115383   3355817975965/11060851509271    0           0       0
                        647845179188/3216320057751      73281519250/8382639484533       552539513391/3454668386233      3354512671639/8306763924573     4040/17871  0       0
                        82889/524892                    0                               15625/83664                     69875/102672                    -2260/8211  0.25    0   ]
        RK.order    = 4
        RK.GSA          = false
        RK.const_diag   = true
    elseif RKMETHOD=="KC4B"
        #4-th order Kennedy-Carpenter ARK4(3)7L[2]SA1 from ANM2019:
        RK.AI, RK.AE= ARK4_3_7L_2_SA1()
        RK.order    = 4
        av          = [ 0.0, 
                        0.0, 
                        4469248916618/8635866897933, 
                        -621260224600/4094290005349,
                        696572312987/2942599194819,
                        1532940081127/5565293938103,
                        2441/20000, 
                        0.0                         ]
        RK.GSA          = false
        RK.const_diag   = true
	elseif RKMETHOD=="RALSTON4"
		RK.AE       = [ 0 		 	0 				0 			0 			0
						0.4 		0.0				0 			0  			0
						0.29697761  0.15875964 		0.0 		0.0  		0.0
						0.21810040 	-3.05096516 	3.83286476 	0.0  		0.0
						0.17476028 	-0.55148066 	1.20553560 	0.17118478  0.0 ]
        RK.AI       = copy(RK.AE)
        RK.order    = 4
    elseif RKMETHOD=="EX4"
		RK.AE       = [ 0 		 	0 				0 			0 			0
						0.5 		0.0				0 			0  			0
						0.0         0.5 		    0.0 		0.0  		0.0
						0.0         0.0             1.0         0.0         0.0 
						1/6 	    1/3	            1/3 	    1/6  		0.0 ]
        RK.AI       = copy(RK.AE)
        RK.order    = 4
    elseif RKMETHOD=="CALVO4"
        RK.AI       = [ 0.0         0.0         0.0         0.0         0.0         0.0         0.0     
                        0.0         0.25        0.0         0.0         0.0         0.0         0.0
                        0.0         0.5         0.25        0.0         0.0         0.0         0.0     
                        0.0         17/50       -1/25       0.25        0.0         0.0         0.0             
                        0.0         371/1360    -137/2720   15/544      0.25        0.0         0.0             
                        0.0         25/24       -49/48      125/16      -85/12      0.25        0.0     
                        0.0         25/24       -49/48      125/16      -85/12      0.25        0.0  ]
        RK.AE       = [ 0.0         0.0         0.0         0.0         0.0         0.0         0.0  
                        0.25        0.0         0.0         0.0         0.0         0.0         0.0 
                        -0.25       1.0         0.0         0.0         0.0         0.0         0.0     
                        -13/100     43/75       8/75        0.0         0.0         0.0         0.0     
                        -6/85       42/85       179/1360    -15/272     0.0         0.0         0.0     
                        0.0         79/24       -5/8        25/2        -85/6       0.0         0.0     
                        0.0         25/24       -49/48      125/16      -85/12      0.25        0.0  ]
        RK.order    = 4
        RK.GSA          = false
        RK.const_diag   = true
    elseif RKMETHOD=="TSZ3"
        RK.AE       = [ 0.0                 0.0                 0.0                 0.0                 0.0     
                        0.520300000000000   0.0                 0.0                 0.0                 0.0 
                        0.026500000000000   0.938000000000000   0.0                 0.0                 0.0 
                        0.122175553766880   0.105600000000000   0.018300000000000   0.0                 0.0 
                        -0.033950868284890  0.218016324016351   0.258600000000000   0.557334544268539   0.0 ]
        RK.AI       = [ 0.0                 0.0                 0.0                 0.0                 0.0
                        -0.520300000000000  0.520300000000000   0.0                 0.0                 0.0
                        0.911500000000000   -1.876000000000000  0.964500000000000   0.0                 0.0
                        -0.401069249711528  0.663393695944647   -0.508400000000000  0.246075553766880   0.0
                        -0.155925222099085  -0.084089256959580  -1.070724285228281  0.310738764286946   1.0 ] + RK.AE
        RK.order    = 3
        RK.GSA          = true
        RK.const_diag   = false
    elseif RKMETHOD=="CLASSICRK4_E"
        RK.AE           = [ 0.0     0.0     0.0     0.0     0.0
                            0.5     0.0     0.0     0.0     0.0
                            0.0     0.5     0.0     0.0     0.0
                            0.0     0.0     1.0     0.0     0.0
                            1/6     1/3     1/3     1/6     0.0 ]
        RK.AI           = copy(RK.AE)
        RK.order        = 4
        RK.GSA          = false
        RK.const_diag   = true
    elseif RKMETHOD=="CLASSICRK3_I"
        RK.AI           = [ 0.0     0.0     0.0     0.0     0.0
                            0.0     1/2     0.0     0.0     0.0
                            0.0     1/6     1/2     0.0     0.0
                            0.0     -1/2    1/2     1/2     0.0
                            0.0     3/2     -3/2    1/2     1/2 ]
        RK.AE           = NaN*RK.AI
        RK.order        = 3
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD=="CLASSICRK4_I"
        RK.AI           = [ 0.0         0.0         0.0         0.0     0.0     0.0
                            0.0         1/4         0.0         0.0     0.0     0.0
                            0.0         1/2         1/4         0.0     0.0     0.0
                            0.0         17/50       -1/25       1/4     0.0     0.0 
                            0.0         371/1360    -137/2720   15/544  1/4     0.0
                            0.0         25/24       -49/48      125/16  -85/12  1/4 ]
        RK.AE           = NaN*RK.AI
        RK.order        = 4
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD=="KC35"
        RK.AI           = ESDIRK3_2_5L_2_SA()
        RK.AE           = NaN*RK.AI
        RK.order        = 3
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD=="KC47"
        RK.AI           = ESDIRK4_3_7L_2_SA()
        RK.AE           = NaN*RK.AI
        RK.order        = 4
        RK.GSA          = true
        RK.const_diag   = true
    elseif RKMETHOD=="KC58"
        #Not I-stable:
        RK.AI           = ESDIRK5_4_8L_2_SA()
        RK.AE           = NaN*RK.AI
        RK.order        = 5
        RK.GSA          = true
        RK.const_diag   = true
	else
		error("Unknown RK method=$(RKMethod)")
    end
    
    #Obtain more info from AI and AE:
    RK.stages       = size(RK.AI,1)
    RK.c            = sum(RK.AI,dims=2)[:]

    return RK
    
end

#Plot stability region for ydot = g(y):
function IRK_PlotStabilityRegion(RKMethod::String, xv::Array{Float64,1}, yv::Array{Float64,1})

    #Load coefficients:
    RK      = RK_Coefficients(RKMethod)

    #Compute stability function:
    ss      = RK.stages
    Im      = eye(ss)
    ev      = ones(ss)
    function Rfun(z::Number)
        return (Im-z*RK.AI)\ev
    end
    
    #Plot stability region:
    PlotStabilityRegion(Rfun, ss, xv, yv)
    
    return
    
end
function IRK_PlotStabilityRegion_LO(RKMethod::String, xv::Array{Float64,1}, yv::Array{Float64,1})

    #Load coefficients:
    RK      = RK_Coefficients(RKMethod)

    #Compute stability function:
    ss      = RK.stages
    Im      = eye(ss)
    ev      = ones(ss)
    function Rfun(z::Number)
        return [1.0+z*(transpose(RK.bhatI)*((Im-z*RK.AI)\ev))]
    end
    
    #Plot stability region:
    PlotStabilityRegion(Rfun, 1, xv, yv)
    
    return
    
end

#Plot stability region for ydot = [g(y)-mu*y] + mu*y:
function EINRK_PlotStabilityRegion!(RKMethod::String, muDeltat::Union{Float64,Complex{Float64}}, xv::Array{Float64,1}, yv::Array{Float64,1};
    color::String="b")

    #Load coefficients:
    RK      = RK_Coefficients(RKMethod)

    #Compute stability function:
    ss      = RK.stages
    Im      = eye(ss)
    ev      = ones(ss)
    Rfun    = (z::Number) -> (Im-muDeltat*RK.AI-(z-muDeltat)*RK.AE)\ev
    
    #Plot stability region:
    PlotStabilityRegion!(Rfun, ss, xv, yv, color=color)
    
    return
    
end
function EINRK_PlotStabilityRegion_LO!(RKMethod::String, muDeltat::Union{Float64,Complex{Float64}}, xv::Array{Float64,1}, yv::Array{Float64,1};
    color::String="b")

    #Load coefficients:
    RK      = RK_Coefficients(RKMethod)

    #Compute stability function:
    ss      = RK.stages
    Im      = eye(ss)
    ev      = ones(ss)
    function Rfun_EINRK_LO(z::Number)
        y_stages    = (Im-muDeltat*RK.AI-(z-muDeltat)*RK.AE)\ev
        y_np1       = 1.0 + muDeltat*dot(RK.bhatI, y_stages) + (z-muDeltat)*dot(RK.bhatE, y_stages)
        return [y_np1]
    end
    
    #Plot stability region:
    PlotStabilityRegion!(Rfun_EINRK_LO, 1, xv, yv, color=color)
    
    return
    
end

function PlotStabilityRegion(Rfun::Function, stages::Int, xv::Array{Float64,1}, yv::Array{Float64,1};
    color::String="b")

    figure()
    PlotStabilityRegion!(Rfun, stages, xv, yv, color=color)
    return
    
end

function PlotStabilityRegion!(Rfun::Function, stages::Int, xv::Array{Float64,1}, yv::Array{Float64,1};
    color::String="b")

    #Compute R at each point of the grid:
    xm,ym   = ndgrid(xv,yv)
    Rm      = zeros(Complex,length(xv),length(yv),stages)
    for jj=1:length(yv), ii=1:length(xv)
        z               = xv[ii]+1im*yv[jj]
        Rm[ii,jj,:]     = Rfun(z)
    end
    
    #Plot stability region for each stage:
#     colorv  = PyPlotColors("jet2",stages)
    leg     = []
    for ss=stages
#         PyPlot.contour(xm,ym,abs.(Rm[:,:,ss]),[1.0],colors=[colorv[ss]])
        PyPlot.contourf(xm,ym,abs.(Rm[:,:,ss]),[0.0, 1.0],colors=[color])
#         PyPlot.contour(xm,ym,abs.(Rm[:,:,ss]),[1.1],colors=[color],linestyle="dashed")
#         leg     = vcat(leg,"s=$ss")
#         surf(xm, ym, log10.(abs.(Rm[:,:,ss])), cmap="jet")
#         xlabel("x")
#         ylabel("y")
    end
    # PyPlot.axis("equal")
#     PyPlot.grid(true)
    # PyPlot.legend(leg,loc="best")

    #=
    #Display R(z->-Inf) for each stage:
    Nz      = 100
    zv      = -logspace(0.01, 1e10, Nz); 
    Rm      = zeros(Complex,Nz,stages)
    for ii=1:Nz
        z           = zv[ii]
        Rm[ii,:]    .= Rfun(z)
    end
    figure()
    for ii=1:stages
        loglog(abs.(zv), abs.(Rm[:,ii]), color=colorv[ii])
    end
    title("Real axis")
    
    #Display R(z->-Inf) for each stage:
    Nz      = 100
    zv      = logspace(0.01, 1e10, Nz); 
    Rm      = zeros(Complex,Nz,stages)
    for ii=1:Nz
        z           = 1im*zv[ii]
        Rm[ii,:]    .= Rfun(z)
    end
    figure()
    for ii=1:stages
        loglog(abs.(zv), abs.(Rm[:,ii]), color=colorv[ii])
    end
    title("Imaginary axis")
#     display()
    =#
    
end
