# joJittshift - generates the sampling operator for time-jittered OBC acquisition, for one source vessel with two airgun arrays
# julia version: Yijun Zhang 2020.08.11

## helper module
module joJittshift_etc
    using JOLI
    # planned
    # forward Jittshift 
    function Jittshiftfwd(x::Array, wn::Array, dim::Array, jitdim::Array, gridsjitb1arr1IND::Array, gridsjitb1arr2IND::Array, tfirejitb1arr1gridIND::Array, tfirejitb1arr2gridIND::Array, tshiftb1arr1::Array, tshiftb1arr2::Array)
        Ft = joDFT(dim[1];DDT=Float64);
    
        x = reshape(x,Tuple(dim));
        y1 = zeros(Tuple(jitdim));
        y2 = copy(y1);
        
        for j = 1 : length(gridsjitb1arr1IND)
            shot_arr1 = x[:, :, Int(gridsjitb1arr1IND[j])];
            shot_arr1 = Ft*shot_arr1;
            for k = 1:size(shot_arr1,2)
                shot_arr1[:,k] = shot_arr1[:,k].*exp.(-1im*wn*tshiftb1arr1[j]);
            end
            shot_arr1 = Ft'*shot_arr1;
    
            shot_arr2 = x[:, :, Int(gridsjitb1arr2IND[j])];
            shot_arr2 = Ft*shot_arr2;
            for k = 1:size(shot_arr2,2)
                shot_arr2[:,k] = shot_arr2[:,k].*exp.(-1im*wn*tshiftb1arr2[j]);
            end
            shot_arr2 = Ft'*shot_arr2;
    
            y1[Int(tfirejitb1arr1gridIND[j]) : Int(tfirejitb1arr1gridIND[j] + dim[1] - 1), :] = y1[Int(tfirejitb1arr1gridIND[j]) : Int(tfirejitb1arr1gridIND[j] + dim[1] - 1), :] + shot_arr1;
            y2[Int(tfirejitb1arr2gridIND[j]) : Int(tfirejitb1arr2gridIND[j] + dim[1] - 1), :] = y2[Int(tfirejitb1arr2gridIND[j]) : Int(tfirejitb1arr2gridIND[j] + dim[1] - 1), :] + shot_arr2;
        end
         
        y = y1 + y2;
        y = y[:]; 
        
        return y
    end
    # adjoint Jittshift
    function Jittshiftadj(x::Array, wn::Array, dim::Array, jitdim::Array, gridsjitb1arr1IND::Array, gridsjitb1arr2IND::Array, tfirejitb1arr1gridIND::Array, tfirejitb1arr2gridIND::Array, tshiftb1arr1::Array, tshiftb1arr2::Array)
        Ft = joDFT(dim[1];DDT=Float64);

        x1 = reshape(x, Tuple(jitdim));
        x2 = copy(x1);
        y1 = zeros(Tuple(dim));
        y2 = copy(y1);
        
        for j = 1 : length(gridsjitb1arr1IND)     
          y1[:,:,Int(gridsjitb1arr1IND[j])] = x1[Int(tfirejitb1arr1gridIND[j]) : Int(tfirejitb1arr1gridIND[j] + dim[1] - 1), :];
          shot_arr1 = Ft*y1[:,:,Int(gridsjitb1arr1IND[j])];
          for k = 1:size(shot_arr1,2)
              shot_arr1[:,k] = shot_arr1[:,k].*exp.(-1im*wn*-tshiftb1arr1[j]);
          end   
          y1[:,:,Int(gridsjitb1arr1IND[j])] = Ft'*shot_arr1; 
          
          y2[:,:,Int(gridsjitb1arr2IND[j])] = x2[Int(tfirejitb1arr2gridIND[j]) : Int(tfirejitb1arr2gridIND[j] + dim[1] - 1), :];
          shot_arr2 = Ft*y2[:,:,Int(gridsjitb1arr2IND[j])];
          for k = 1:size(shot_arr2,2)
              shot_arr2[:,k] = shot_arr2[:,k].*exp.(-1im*wn*-tshiftb1arr2[j]);
          end   
          y2[:,:,Int(gridsjitb1arr2IND[j])] = Ft'*shot_arr2; 
        end 
    
        y = y1 + y2;
        y = y[:];

        return y
    end
    # planned =  false
    # forward Jittshift 
    function Jittshiftfwd(x::Array, wn::Array, dim::Array, jitdim::Array, gridsjitb1arr1IND::Array, gridsjitb1arr2IND::Array, tfirejitb1arr1gridIND::Array, tfirejitb1arr2gridIND::Array, tshiftb1arr1::Array, tshiftb1arr2::Array, planned)
        nvc = size(x,2);
        Ft = joDFT(dim[1];DDT=Float64);
    
        x = reshape(x,(dim[1],dim[3],nvc));
        x = permutedims(x,[1,3,2]);
        
        y1 = zeros(jitdim[1],nvc);
        y2 = copy(y1);
        
        for j = 1 : length(gridsjitb1arr1IND)
            shot_arr1 = x[:, :, Int(gridsjitb1arr1IND[j])];
            shot_arr1 = Ft*shot_arr1;
            for k = 1:size(shot_arr1,2)
                shot_arr1[:,k] = shot_arr1[:,k].*exp.(-1im*wn*tshiftb1arr1[j]);
            end
            shot_arr1 = Ft'*shot_arr1;
    
            shot_arr2 = x[:, :, Int(gridsjitb1arr2IND[j])];
            shot_arr2 = Ft*shot_arr2;
            for k = 1:size(shot_arr2,2)
                shot_arr2[:,k] = shot_arr2[:,k].*exp.(-1im*wn*tshiftb1arr2[j]);
            end
            shot_arr2 = Ft'*shot_arr2;
    
            y1[Int(tfirejitb1arr1gridIND[j]) : Int(tfirejitb1arr1gridIND[j] + dim[1] - 1), :] = y1[Int(tfirejitb1arr1gridIND[j]) : Int(tfirejitb1arr1gridIND[j] + dim[1] - 1), :] + shot_arr1;
            y2[Int(tfirejitb1arr2gridIND[j]) : Int(tfirejitb1arr2gridIND[j] + dim[1] - 1), :] = y2[Int(tfirejitb1arr2gridIND[j]) : Int(tfirejitb1arr2gridIND[j] + dim[1] - 1), :] + shot_arr2;
        end
         
        y = y1 + y2;
        y = y[:]; 
        
        return y
    end
    # adjoint Jittshift
    function Jittshiftadj(x::Array, wn::Array, dim::Array, jitdim::Array, gridsjitb1arr1IND::Array, gridsjitb1arr2IND::Array, tfirejitb1arr1gridIND::Array, tfirejitb1arr2gridIND::Array, tshiftb1arr1::Array, tshiftb1arr2::Array, planned)
        nvc = size(x,2)
        Ft = joDFT(dim[1];DDT=Float64);

        #x1 = reshape(x, Tuple(jitdim));
        x1 = x;
        x2 = copy(x1);
        y1 = zeros(dim[1],nvc,dim[3]);
        y2 = copy(y1);
        
        for j = 1 : length(gridsjitb1arr1IND)     
          y1[:,:,Int(gridsjitb1arr1IND[j])] = x1[Int(tfirejitb1arr1gridIND[j]) : Int(tfirejitb1arr1gridIND[j] + dim[1] - 1), :];
          shot_arr1 = Ft*y1[:,:,Int(gridsjitb1arr1IND[j])];
          for k = 1:size(shot_arr1,2)
              shot_arr1[:,k] = shot_arr1[:,k].*exp.(-1im*wn*-tshiftb1arr1[j]);
          end   
          y1[:,:,Int(gridsjitb1arr1IND[j])] = Ft'*shot_arr1; 
          
          y2[:,:,Int(gridsjitb1arr2IND[j])] = x2[Int(tfirejitb1arr2gridIND[j]) : Int(tfirejitb1arr2gridIND[j] + dim[1] - 1), :];
          shot_arr2 = Ft*y2[:,:,Int(gridsjitb1arr2IND[j])];
          for k = 1:size(shot_arr2,2)
              shot_arr2[:,k] = shot_arr2[:,k].*exp.(-1im*wn*-tshiftb1arr2[j]);
          end   
          y2[:,:,Int(gridsjitb1arr2IND[j])] = Ft'*shot_arr2; 
        end 
    
        y = y1 + y2;
        y = y[:];

        return y
    end
end
using .joJittshift_etc

export joJittshift
"""
    julia> joJittshift(dim,jitdim,jitacq,wn)
    joJittshift generates the sampling operator for time-jittered OBC acquisition, for one source vessel with two airgun arrays
# Signature
    joJittshift(dim::Int, jitdim::Int, jitacq::Dict, wn::Array;DDT::DataType=joFloat,RDT::DataType=joFloat)
# Arguments
- dim : ?
- jitdim : ?
- jitacq : ?
- wn : ?
# Notes
- ?
# Examples
- ?
"""
function joJittshift(dim::Array, jitdim::Array, jitacq::Dict, wn::Array;planned::Bool=true, DDT::DataType=joFloat,RDT::DataType=joFloat)
    
    # Jittered acquisition parameters
    gridsjitb1arr1IND = jitacq["gridsjitb1arr1IND"];
    gridsjitb1arr2IND = jitacq["gridsjitb1arr2IND"];
    tfirejitb1arr1gridIND = jitacq["tfirejitb1arr1gridIND"];
    tfirejitb1arr2gridIND = jitacq["tfirejitb1arr2gridIND"];
    tshiftb1arr1 = jitacq["tshiftb1arr1"];
    tshiftb1arr2 = jitacq["tshiftb1arr2"];

    if planned 
        # Size of the sampling operator
        m = prod(jitdim);
        n = prod(dim);

        joLinearFunctionFwd_A(m,n,
            v1->joJittshift_etc.Jittshiftfwd(v1, wn, dim, jitdim, gridsjitb1arr1IND, gridsjitb1arr2IND, tfirejitb1arr1gridIND, tfirejitb1arr2gridIND, tshiftb1arr1, tshiftb1arr2),
            v2->joJittshift_etc.Jittshiftadj(v2, wn, dim, jitdim, gridsjitb1arr1IND, gridsjitb1arr2IND, tfirejitb1arr1gridIND, tfirejitb1arr2gridIND, tshiftb1arr1, tshiftb1arr2),
            DDT,RDT;
            name="joJittshift")
    else #planned == false, have to parallel the data
        # Size of the sampling operator
        m = jitdim[1];
        n = dim[1]*dim[3];

        joLinearFunctionFwd_A(m,n,
            v1->joJittshift_etc.Jittshiftfwd(v1, wn, dim, jitdim, gridsjitb1arr1IND, gridsjitb1arr2IND, tfirejitb1arr1gridIND, tfirejitb1arr2gridIND, tshiftb1arr1, tshiftb1arr2, planned),
            v2->joJittshift_etc.Jittshiftadj(v2, wn, dim, jitdim, gridsjitb1arr1IND, gridsjitb1arr2IND, tfirejitb1arr1gridIND, tfirejitb1arr2gridIND, tshiftb1arr1, tshiftb1arr2, planned),
            DDT,RDT;
            name="joJittshift")
    end
end

