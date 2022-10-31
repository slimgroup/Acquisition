export Time_Jitt_LR

function Time_Jitt_LR(A, x,g,params)

    # julia version: Yijun Zhang 2020.09.15


    # Use:
    #   Time_Jitt_LR(x,g,params)
    #
    # Input:   
    #      x - unknow data 
    #      g - computed gradient       
    #      params - parameter file
    
    # Author: Rajiv Kumar
    #         Seismic Laboratory for Imaging and Modeling
    #         Department of Earth, Ocean, and Atmospheric Sciences
    #         The University of British Columbia
    #         
    # Date: January, 2015
    
    # You may use this code only under the conditions and terms of the
    # license contained in the file LICENSE provided with this source
    # code. If you do not agree to these terms you may not use this
    # software.
    #----------------------------------------------------------------------------------------------------
    typeof(params["nt"]) == Int ? params["nt"] = params["nt"] : params["nt"] = Int(params["nt"])
    typeof(params["nr"]) == Int ? params["nr"] = params["nr"] : params["nr"] = Int(params["nr"])
    typeof(params["nc"]) == Int ? params["nc"] = params["nc"] : params["nc"] = Int(params["nc"])
    typeof(params["nf"]) == Int ? params["nf"] = params["nf"] : params["nf"] = Int(params["nf"])
    typeof(params["nm"]) == Int ? params["nm"] = params["nm"] : params["nm"] = Int(params["nm"])
    typeof(params["nh"]) == Int ? params["nh"] = params["nh"] : params["nh"] = Int(params["nh"])
    typeof(params["k"]) == Int ? params["k"] = params["k"] : params["k"] = Int(params["k"])
    
    # replacement of NLfunforward (this define complete A)
    Ft = joDFTR(params["nt"];DDT = Float32);
    #Ft = joKron(joDirac(params["nr"]*params["nc"]; DDT= Complex{Float64}),Ft); # this operator need to be parallized 
    L = x[1:params["nf"]*params["nm"]*params["k"]];
    R = x[params["nf"]*params["nm"]*params["k"]+1:end];
    L = convert(Array{Complex{Float32},3},reshape(L,params["nf"],params["nm"],params["k"]));
    R = convert(Array{Complex{Float32},3},reshape(R,params["nf"],params["nh"],params["k"]));
    #using SeisJOLI
    MH = joSRtoCMO(params["nr"],params["nc"];DDT = Complex{Float32});  
    
    
    if isempty(g)
        output = zeros(Complex{Float32},params["nf"],params["nr"]*params["nc"]);  #this one is the distributed array
        for i = 1:params["nf"]
            output[i,:] = MH'*vec(L[i,:,:]*R[i,:,:]'); # set array to Darray (remotecall)
        end
        
        f1 = params["RM"]*Float64.(vec(Ft'*output)); output = nothing; #Ft and Rm should be Distributed operator
        f2 = 0;
    else
        fp  = Ft*Float32.(reshape(params["RM"]'*vec(g),params["nt"],params["nr"]*params["nc"])); g = nothing;# g should be Darray
        fp  = reshape(fp,params["nf"],params["nr"],params["nc"]); #fp is the Darray
        
        fL  = zeros(Complex{Float32},params["nf"],params["nm"],params["k"]);
        fR  = zeros(Complex{Float32},params["nf"],params["nh"],params["k"]);
        fLR = zeros(Complex{Float32},params["nf"],params["nm"],params["nh"]); # this one need to be parallized 
        
        for i = 1:params["nf"]
            ftest = reshape(MH*vec(fp[i,:,:]),params["nm"],params["nh"]);
            fLR[i,:,:] = ftest;
            fL[i,:,:]  = ftest*R[i,:,:]; #get local Darray and set Darray to Array
            fR[i,:,:]  = ftest'*L[i,:,:];
        end
	
        f1 = [vec(fL); vec(fR)];
        f2 = vec(fLR);
    end
    return f1, f2
end
    
