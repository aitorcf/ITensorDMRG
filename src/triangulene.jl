using ITensors
using LinearAlgebra
using DelimitedFiles

triangulenedir = joinpath(dirname(pathof(ITensorDMRG)), "..", "triangulene" )

qn_vector = [
 QN(("Nz",-1,-1),("Sz",-3)) => 4,
 QN(("Nz",-1,-1),("Sz",-1)) => 24,
  QN(("Nz",-1,-1),("Sz",1)) => 24,
  QN(("Nz",-1,-1),("Sz",3)) => 4,
  QN(("Nz",0,-1),("Sz",-2)) => 16,
   QN(("Nz",0,-1),("Sz",0)) => 36,
   QN(("Nz",0,-1),("Sz",2)) => 16,
  QN(("Nz",1,-1),("Sz",-3)) => 4,
  QN(("Nz",1,-1),("Sz",-1)) => 24,
   QN(("Nz",1,-1),("Sz",1)) => 24,
   QN(("Nz",1,-1),("Sz",3)) => 4
]
function ITensors.space(::SiteType"Triangulene";conserve_qns=false)
    if conserve_qns 
        return qn_vector
    else
        return 180
    end
end

sz_matrix = diagm(reduce( vcat , Float64[val(qn,"Sz")/2.0 for _ in 1:dim] for (qn,dim) in qn_vector ))
ITensors.op(::OpName"Sz",::SiteType"Triangulene") = sz_matrix

spectrum_matrix = diagm(readdlm("$(triangulenedir)/States.txt",Float64)[:,3])
ITensors.op(::OpName"Spectrum",::SiteType"Triangulene") = spectrum_matrix

cdagup1_matrix = readdlm("$(triangulenedir)/Creation_op_up_1.txt",Float64)
cup1_matrix = adjoint(cdagup1_matrix)
ITensors.op(::OpName"Cdagup1",::SiteType"Triangulene") = cdagup1_matrix
ITensors.op(::OpName"Cup1",::SiteType"Triangulene") = cup1_matrix

cdagdn1_matrix = readdlm("$(triangulenedir)/Creation_op_down_1.txt",Float64)
cdn1_matrix = adjoint(cdagdn1_matrix)
ITensors.op(::OpName"Cdagdn1",::SiteType"Triangulene") = cdagdn1_matrix
ITensors.op(::OpName"Cdn1",::SiteType"Triangulene") = cdn1_matrix

cdagup2_matrix = readdlm("$(triangulenedir)/Creation_op_up_2.txt",Float64)
cup2_matrix = adjoint(cdagup2_matrix)
ITensors.op(::OpName"Cdagup2",::SiteType"Triangulene") = cdagup2_matrix
ITensors.op(::OpName"Cup2",::SiteType"Triangulene") = cup2_matrix

cdagdn2_matrix = readdlm("$(triangulenedir)/Creation_op_down_2.txt",Float64)
cdn2_matrix = adjoint(cdagdn2_matrix)
ITensors.op(::OpName"Cdagdn2",::SiteType"Triangulene") = cdagdn2_matrix
ITensors.op(::OpName"Cdn2",::SiteType"Triangulene") = cdn2_matrix

cdagup3_matrix = readdlm("$(triangulenedir)/Creation_op_up_3.txt",Float64)
cup3_matrix = adjoint(cdagup3_matrix)
ITensors.op(::OpName"Cdagup3",::SiteType"Triangulene") = cdagup3_matrix
ITensors.op(::OpName"Cup3",::SiteType"Triangulene") = cup3_matrix

cdagdn3_matrix = readdlm("$(triangulenedir)/Creation_op_down_3.txt",Float64)
cdn3_matrix = adjoint(cdagdn3_matrix)
ITensors.op(::OpName"Cdagdn3",::SiteType"Triangulene") = cdagdn3_matrix
ITensors.op(::OpName"Cdn3",::SiteType"Triangulene") = cdn3_matrix

cdagup4_matrix = readdlm("$(triangulenedir)/Creation_op_up_4.txt",Float64)
cup4_matrix = adjoint(cdagup4_matrix)
ITensors.op(::OpName"Cdagup4",::SiteType"Triangulene") = cdagup4_matrix
ITensors.op(::OpName"Cup4",::SiteType"Triangulene") = cup4_matrix

cdagdn4_matrix = readdlm("$(triangulenedir)/Creation_op_down_4.txt",Float64)
cdn4_matrix = adjoint(cdagdn4_matrix)
ITensors.op(::OpName"Cdagdn4",::SiteType"Triangulene") = cdagdn4_matrix
ITensors.op(::OpName"Cdn4",::SiteType"Triangulene") = cdn4_matrix

nup1_matrix = cdagup1_matrix*cup1_matrix
ndn1_matrix = cdagdn1_matrix*cdn1_matrix
n1_matrix = nup1_matrix .+ ndn1_matrix
ITensors.op(::OpName"N1up",::SiteType"Triangulene") = nup1_matrix
ITensors.op(::OpName"N1dn",::SiteType"Triangulene") = ndn1_matrix
ITensors.op(::OpName"N1",::SiteType"Triangulene") = n1_matrix

nup2_matrix = cdagup2_matrix*cup2_matrix
ndn2_matrix = cdagdn2_matrix*cdn2_matrix
n2_matrix = nup2_matrix .+ ndn2_matrix
ITensors.op(::OpName"N2up",::SiteType"Triangulene") = nup2_matrix
ITensors.op(::OpName"N2dn",::SiteType"Triangulene") = ndn2_matrix
ITensors.op(::OpName"N2",::SiteType"Triangulene") = n2_matrix

nup3_matrix = cdagup3_matrix*cup3_matrix
ndn3_matrix = cdagdn3_matrix*cdn3_matrix
n3_matrix = nup3_matrix .+ ndn3_matrix
ITensors.op(::OpName"N3up",::SiteType"Triangulene") = nup3_matrix
ITensors.op(::OpName"N3dn",::SiteType"Triangulene") = ndn3_matrix
ITensors.op(::OpName"N3",::SiteType"Triangulene") = n3_matrix

nup4_matrix = cdagup4_matrix*cup4_matrix
ndn4_matrix = cdagdn4_matrix*cdn4_matrix
n4_matrix = nup4_matrix .+ ndn4_matrix
ITensors.op(::OpName"N4up",::SiteType"Triangulene") = nup4_matrix
ITensors.op(::OpName"N4dn",::SiteType"Triangulene") = ndn4_matrix
ITensors.op(::OpName"N4",::SiteType"Triangulene") = n4_matrix

n_matrix = n1_matrix .+ n2_matrix .+ n3_matrix .+ n4_matrix
ITensors.op(::OpName"N",::SiteType"Triangulene") = n_matrix
