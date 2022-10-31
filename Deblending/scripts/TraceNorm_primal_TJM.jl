export TraceNorm_primal_TJM

function TraceNorm_primal_TJM(x, weights, params)
    # Input: x-unknown data; weight-weight value; params-parameter file
    p = 0.5 * norm(x.*weights)^2;
    return p
end
