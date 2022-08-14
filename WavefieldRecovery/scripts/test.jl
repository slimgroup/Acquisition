using DrWatson
@quickactivate "WavefieldRecovery"
using MAT
include(projectdir()*"/scripts/SubFunctions/ConventionalLR.jl")
include(projectdir()*"/scripts/SubFunctions/WeightLR.jl")
#load the true data
datafile  = MAT.matopen(projectdir()*"/data/Full.mat")
data_full = read(datafile,"Full_Data")
datafull = copy(data_full[:,1:355,1:355])
close(datafile)
nff,nr,ns = size(datafull)

#load the jittered indexed for subsampling
sub_index = MAT.matopen(projectdir()*"/data/Ind.mat")
ind = read(sub_index, "ind")
close(sub_index)

###define operators subsampling and joMH(source-receiver => midpoint-offset)
RM = joKron(joRestriction(ns,dropdims(trunc.(Int,ind),dims=1);DDT = Complex{Float64}), joDirac(nr;DDT = Complex{Float64}))
MH = joSRtoCMO(nr,ns; DDT = Complex{Float64})

###define rank information if rank=r_sub(convertional weighted method)
rank1 = 30;

###define the final result to save
Final_result = zeros(Complex{Float64},nff,nr,ns)
pre_L = [];
pre_R = [];
###conventional LR reconstruction
for SliceNum = 123:124 
    
    if SliceNum == 123
        ### created subsampling data
        SubD = reshape(RM'*(RM*vec(data_full[SliceNum,:,:])),nr,ns)
        ### reconstructed data via conventional LR
        L, R = ConventionalLR(SubD, MH, rank1; Df_test = data_full[SliceNum,:,:]);
        Final_result[SliceNum, :, :] = reshape(MH'*vec(L*R'), nr,ns)
        global pre_L = copy(L);
        global pre_R = copy(R);
    else 
        ### created subsampling data
        SubD = reshape(RM'*(RM*vec(data_full[SliceNum,:,:])),nr,ns)
        ### reconstructed data via weighted LR
        L, R = WeightedLR(pre_L, pre_R, SubD, MH, rank1; Df_test = data_full[SliceNum,:,:])
        Final_result[SliceNum, :, :] = reshape(MH'*vec(L*R'), nr,ns)
    end
end



