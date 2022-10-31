export TraceNorm_dual_TJM

function TraceNorm_dual_TJM(x, weights, params)
    # Input: x-unknown data; weight-weight value; params-parameter file
    # dual of trace norm is operator norm i.e maximum singular value

    typeof(params["nf"]) == Int ? params["nf"] = params["nf"] : params["nf"] = Int(params["nf"])
    typeof(params["nm"]) == Int ? params["nm"] = params["nm"] : params["nm"] = Int(params["nm"])
    typeof(params["nh"]) == Int ? params["nh"] = params["nh"] : params["nh"] = Int(params["nh"])

    x = reshape(x, params["nf"], params["nm"],params["nh"]);
    dfull = zeros(params["nf"],1);
    for i in 1:params["nf"]
	    tmp = svds(x[i,:,:]; nsv = 1, ritzvec = false)[1]; #using Arpack
        dfull[i] = tmp.S[1];
    end
    d = maximum(dfull);
    dfull = nothing
    return d
end
