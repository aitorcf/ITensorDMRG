using Pkg; Pkg.activate("/home/aitor/Bulegoa/ITensorDMRG")

using ITensorDMRG
using ITensors
using LinearAlgebra

let

    N = 6
    sites = siteinds("Triangulene",N;conserve_qns=true)

    h = 10.0
    mu = -4.0

    # all together
    os = OpSum()

    # on-site part
    os_onsite = OpSum()
    for i in 1:N

        os_onsite += "Spectrum",i
        os_onsite += mu,"N",i

    end
    H_onsite = MPO(os_onsite,sites)

    # hopping part
    os_hop = OpSum()
    for i in 1:(N-1)

        os_hop += h,"Cdagup1",i,"Cup1",i+1
        os_hop += h,"Cup1",i,"Cdagup1",i+1
        os_hop += h,"Cdagdn1",i,"Cdn1",i+1
        os_hop += h,"Cdn1",i,"Cdagdn1",i+1

        os_hop += h,"Cdagup2",i,"Cup2",i+1
        os_hop += h,"Cup2",i,"Cdagup2",i+1
        os_hop += h,"Cdagdn2",i,"Cdn2",i+1
        os_hop += h,"Cdn2",i,"Cdagdn2",i+1

        os_hop += h,"Cdagup3",i,"Cup3",i+1
        os_hop += h,"Cup3",i,"Cdagup3",i+1
        os_hop += h,"Cdagdn3",i,"Cdn3",i+1
        os_hop += h,"Cdn3",i,"Cdagdn3",i+1

        os_hop += h,"Cdagup4",i,"Cup4",i+1
        os_hop += h,"Cup4",i,"Cdagup4",i+1
        os_hop += h,"Cdagdn4",i,"Cdn4",i+1
        os_hop += h,"Cdn4",i,"Cdagdn4",i+1

    end
    H_hop = MPO(os_hop,sites)

    # dmrg parameters
    nsweeps = 7
    maxdim = [10,20,20]#,100,200,500,500] # gradually increase states kept
    cutoff = [1E-9] # desired truncation error
    noise = [1E-6, 1E-7, 1E-8, 0.0]

    # initial guess
    guess= reduce(vcat,[[5,29,6,30] for _ in 1:(N÷4+1)])[1:N]
    guess = [73 for _ in 1:N]
    psi0 = randomMPS(sites,guess)

    #psi0 = MPS(sites,guess)
    println( "Guess state:" )
    @show maxlinkdim(psi0)
    println()

    # dmrg iterations
    energy,psi = dmrg([H_hop,H_onsite],psi0; nsweeps, maxdim, cutoff,noise)

    # spin correlations
    corr_Sz = 1.0*correlation_matrix(psi,"Sz","Sz")#;sites=[1,N_sites])
    exp_Sz = expect(psi,"Sz");
    @show flux(psi)
    @show exp_Sz
    @show corr_Sz

end
