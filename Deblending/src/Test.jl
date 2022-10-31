using DrWatson
@quickactivate "Deblending"

# load groundtruth data
using MAT, Random, JOLI, SeisJOLI, GenSPGL, Arpack, LinearAlgebra

include(projectdir()*"/scripts/SubFunctions/jitter_airgunarrays.jl")
include(projectdir()*"/scripts/SubFunctions/joJittshift.jl")
include(projectdir()*"/scripts/SubFunctions/TraceNorm_project_TJM.jl")
include(projectdir()*"/scripts/SubFunctions/TraceNorm_primal_TJM.jl")
include(projectdir()*"/scripts/SubFunctions/TraceNorm_dual_TJM.jl")
include(projectdir()*"/scripts/SubFunctions/joDFTR.jl")
datafile  = MAT.matopen(projectdir()*"/data/D.mat")

data = read(datafile,"D");
close(file)

# number of time, receiver and source
nt,nr,ns = size(data);

# (dt, dr, ds): time, receiver, and source sampling intervals
dt = 0.004;  
dr = 25.0;     
ds = 25.0; 

#Frequency samples
wn = collect([0:nt/2;-nt/2+1:-1]*2*pi/(nt*dt));

# Parameters for time-jittered acquisition
rndfactor = [1000 100];
p = 2;
nboats = 1;
rseed = [8 3];
boatspeed = 2.5;
tfireint_min = 10.0;
tdelay = 10.0;
delayboat = 0;

# Setup the time-jittered (or blended) acquisition
jitacq = jitter_airgunarrays(ns, ds, dt, rndfactor, p, nboats, rseed, boatspeed, tfireint_min, tdelay, delayboat);

# Dimensions for time-jittered acquisition
jitdim = [(length(jitacq["tfirejitgrid"]) + nt - 1)  nr];

# Sampling operator

RM = joJittshift(Int.([nt nr ns]), Int.([jitdim[1] jitdim[2]]), jitacq, wn);

# Generate time-jittered (or blended) data volume [jitD]
JitD = reshape(RM*data[:],Tuple(jitdim));

# observed data
b = JitD[:];

# spgl1 parameter
iteration = 200;
opts = spgOptions( optTol = Float32(1e-5),
                   bpTol = Float32(1e-5),
                   decTol = Float32(1e-6),
                   project = TraceNorm_project_TJM, 
                   primal_norm = TraceNorm_primal_TJM, 
                   dual_norm = TraceNorm_dual_TJM, 
                   proxy = true, 
                   ignorePErr = true, 
                   iterations = iteration,
                   #fid =fid,
                   weights = [1])
include(projectdir()*"/scripts/SubFunctions/Time_Jitt_LR.jl")

rank1 = 50
params = Dict("nm" => ns,
              "nh" => ns*2-1,
              "nf" => floor(Int, nt/2)+1,
              "nr" => nr,
              "nc" => ns,
              "nt" => nt,
              "mode" => 1,
              "RM" => RM,
              "funForward"=> Time_Jitt_LR,
              "ls"=> 1,
              "k"=> rank1)

sigma       = norm(b,2);
LInit       = randn(params["nf"],params["nm"],params["k"])+im*randn(params["nf"],params["nm"],params["k"]);
RInit       = randn(params["nf"],params["nh"],params["k"])+im*randn(params["nf"],params["nh"],params["k"]);
initiweight = 1e-6;
xinit       = initiweight*[vec(LInit);vec(RInit)];
tau         = norm(xinit,1);
sigmaerror = 1e-3;
sigmafact   = sigmaerror*sigma;
# deblending data
xest = spgl1(Time_Jitt_LR,b,x = vec(xinit), tau = tau,sigma = sigmafact,options =opts,params = params);

L    = xest[1][1:params["nf"]*params["nm"]*params["k"]];
R    = xest[1][params["nf"]*params["nm"]*params["k"]+1:end];
L    = reshape(L,params["nf"],params["nm"],params["k"]);
R    = reshape(R,params["nf"],params["nh"],params["k"]);

Ft = joDFTR(params["nt"]);
MH = joSRtoCMO(params["nr"],params["nc"]);  
# reconstruct output
output = zeros(Complex{Float64},params["nf"],params["nr"]*params["nc"]); 
for i = 1:params["nf"]
    output[i,:] = MH'*vec(L[i,:,:]*R[i,:,:]'); 
end

Result = reshape(Ft'*output,nt,nr,ns);

SNR = -20*log10(norm(vec(data)-vec(Result))/norm(vec(data)))
println(SNR)


