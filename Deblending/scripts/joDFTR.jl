module joDFTR_etc
    using FFTW
    using JOLI: jo_convert
    # forward Jittshift
    function DFTRfwd(x::Array, m, n, rdt::DataType)
        nc = size(x,2);
        tmp = fft(x)/sqrt(n);
        if mod(n,2) == 1
            output = tmp[1:m,:]+ [zeros(1,nc); conj(tmp[end:-1:m+1,:])];
        else
            output = tmp[1:m,:] + [zeros(1,nc); conj(tmp[end:-1:m+1,:]); zeros(1,nc)];
        end
        output=jo_convert(rdt,output,false)
        return output
    end
    # adjoint Jittshift
    function DFTRadj(x::Array, m, n, rdt::DataType)
        nc = size(x,2);
        if mod(n,2) == 1
            tmp = [x; conj(x[end:-1:2,:])];
        else
            tmp = [x ; conj(x[end-1:-1:2,:])];
        end
        output = ifft(tmp)*sqrt(n);
        output=jo_convert(rdt,output,false)
        return output
    end
end
using .joDFTR_etc

export joDFTR

function joDFTR(n::Int; DDT::DataType=joFloat, RDT::DataType=(DDT<:Real ? Complex{DDT} : DDT))
    m = Int(floor(n/2) + 1);
    joLinearFunctionFwd_A(m,n,
            v1->joDFTR_etc.DFTRfwd(v1, m, n, RDT),
            v2->joDFTR_etc.DFTRadj(v2, m, n, DDT),
            DDT,RDT;
            name="joDFTR")

end
