export ConventionalLR
using JOLI, GenSPGL, SeisJOLI, Arpack, LinearAlgebra, Random
rng = MersenneTwister(12345);
function ConventionalLR(subD::Tx, MH::joLinearFunction, rank::Int; Df_test::Tx = zeros(ComplexF64, 100,100), r_sub::Int =1)where {ETx<:Number, Tx<:AbstractMatrix{ETx}}
    #subD: Observed data
    #MH: source-receiver <=> midpoint - offset
    #rank: rank information
    #Df_test: ground truth data which is used to subsampling and calculate the signal noise ratio 
    #r_sub: rank for limited subspaces 
    
    r_sub = r_sub==1 ? rank : r_sub;
    nr,ns = size(subD);
        
    ##transfer subsampling data to midoff domain
    midD = reshape(MH*subD[:], nr, 2*ns-1)
    b = copy(midD)
        
    ###prior information
    LInit   = randn(rng, ComplexF32, nr, rank);
    RInit   = randn(rng, ComplexF32, 2*ns-1, rank);
    initfact = 1e-3;
    xinit   = initfact*[vec(copy(LInit));vec(copy(RInit))]

    ###set relative parameters
    tau = norm(copy(xinit[:]),1)
    sigmafact = 1e-2
    sigma   = sigmafact*norm(copy(b[:]),2);

    ###define the parameters for params
    afunT(x) = reshape(x[:],nr,2*ns-1)
    params = Dict("nr"=>rank,  # define rank value here,
                  "nr1"=>r_sub,
                  "Ind"=> vec(midD) .== 0,
                  "numr"=> nr,
                  "numc"=> 2*ns-1,
                  "funForward"=> NLfunForward,
                  "afunT"=> afunT,
                  "afun"=> afun,
                  "mode"=> 1,
                  "ls"=> 1,
                  "logical"=> 0,
                  "funPenalty"=> funLS)
        
    ### Choose the parameters for GenSPGL
    opts = spgOptions(optTol = 1e-5,
                      bpTol = 1e-5,
                      decTol = 1e-4,
                      project = TraceNorm_project,
                      primal_norm = TraceNorm_primal,
                      dual_norm = TraceNorm_dual,
                      proxy = true,
                      ignorePErr = true,
                      iterations = 150,
                      verbosity = 1,
                      funCompositeR = GenSPGL.funCompR1)

    ### using spgl1 to implenment the algorithm
    xLS_jl, r, g, info = spgl1(NLfunForward, b[:], x = xinit[:], tau = tau,
                                sigma = sigma, options = opts, params = params)

    ### extracted the results
    L_Num = params["numr"]*params["nr"]
    L1 = xLS_jl[1:L_Num]
    R1 = xLS_jl[L_Num+1:end]
    L = reshape(L1,params["numr"],params["nr"]);
    R = reshape(R1,params["numc"],params["nr"]);

    ### obtain the final results
    Result = reshape(MH'*vec(L*R'), nr,ns)

    ### calculate the SNR if exsit groundtruth
    if norm(Df_test) != 0
        SNR = -20*log10(norm(vec(Df_test)-vec(Result))/norm(vec(Df_test)));

        println("SNR = ",SNR)
    end
    return L,R

end

