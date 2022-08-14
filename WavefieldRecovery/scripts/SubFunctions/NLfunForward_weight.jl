export NLfunForward_weight


"""
weighted forward function
"""
function NLfunForward_weight(A, x, g, params)
    
    e = params["numr"]*params["nr"]
    w = params["w"]

    e1 = params["numr"]*params["nr1"]
    W_data = params["weight"]
    U_t = copy(W_data[1:e1,1])
    V_t = copy(W_data[e1+1:end,1])
    U = reshape(copy(U_t),params["numr"], params["nr1"])
    V = reshape(copy(V_t),params["numc"], params["nr1"])

    L = x[1:e]
    R = x[(e+1):end]
    L = reshape(L, params["numr"], params["nr"])
    R = reshape(R, params["numc"], params["nr"])

    if isempty(g)
        L += (1/w-1) * U * (U' * L)
        R += (1/w-1) * V * (V' * R);
        f1 = params["afun"](L*(R)', params)
        f2 = 0.
    else
        L += (1/w-1) * U * (U' * L)
        R += (1/w-1) * V * (V' * R);
        fp = params["afunT"](g)
        f1 = [vec((fp + (1/w-1)*U*(U'*fp))*R); vec((fp' + (1/w-1)*V*(V'*fp'))*L)];
        f2_med = fp + (1/w-1)*U*(U'*fp);
        f2 = vec(f2_med+(1/w-1)*f2_med*V*V')
    end

    return f1,f2
end
