using DelimitedFiles
using ITensors

function get_space( spectrumfile::String )

    data = readdlm(spectrumfile)

    qns = Set(data[i,2:3] for i in axes(data,1))
    qndims = Dict(
        qn=>sum( 1 for i in axes(data,1) if data[i,2:3]==qn )
        for qn in qns
    )

    return sort([
        QN(("Nz",Int64(qns[1]),-1),("Sz",Int64(qns[2])))=>dim for (qns,dim) in qndims ],
        by=x->(val(x[1],"Nz"),val(x[1],"Sz"))
    )
end

function get_spectrum( spectrumfile::String )

    data = readdlm(spectrumfile)

    return data[:,1]
end
