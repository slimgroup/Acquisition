export TraceNorm_project_TJM

function TraceNorm_project_TJM(x, B,weights, params)
    # Input: x-unknown data; weight-weight value; params-parameter file; B-tau value
    c = sqrt(B/(0.5*norm(x)^2));			
    x = min(1,c).*x

    #Dummy
    itn = 1	
    return x, itn
end
