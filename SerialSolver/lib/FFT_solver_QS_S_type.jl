include("build_Yle_S.jl")
include("compute_Z_self.jl")
using SparseArrays, SuperLU, IterativeSolvers, FFTW, FourierTools
using LinearMaps, MAT, LinearAlgebra


function FFT_solver_QS_S_type(freq, escalings, incidence_selection, FFTCP, FFTCLp, diagonals, ports, lumped_elements, expansions, GMRES_settings, Zs_info, QS_Rcc_FW)
    # file = matopen("/tmp/matfile.mat", "w")
    # write(file, "frequencies", collect(freq))
    # close(file)
    # FFTCP[2,1] = FFTCP[1,1]
    # FFTCP[3,1] = FFTCP[1,1]
    # FFTCP[3,2] = FFTCP[1,1]
    # FFTCLp[1,2] = FFTCLp[1,1]
    # FFTCLp[2,2] = FFTCLp[1,1]
    # FFTCLp[3,2] = FFTCLp[1,1]
    # lumped_elements.type = convert(Array{Int64}, lumped_elements.type)
    # lumped_elements.value = convert(Array{Float64}, lumped_elements.value)
    # matwrite("/Users/edgardovittoria/Downloads/matfile.mat", Dict(
    #     "freq" => freq,
    #     "escalings" => escalings,
    #     "diagonals" => diagonals,
    #     "ports" => ports,
    #     "lumped_elements" => lumped_elements,
    #     "incidence_selection" => incidence_selection,
    #     "FFTCP" => FFTCP,
    #     "FFTCLp" => FFTCLp,
    #     "expansions" => expansions,
    #     "GMRES_settings" => GMRES_settings,
    #     "Zs_info" => Zs_info,
    #     "QS_Rcc_FW" => QS_Rcc_FW
    # ))
    # vars = matread("/tmp/matfile.mat")
    # println(vars)
    freq = freq .* escalings["freq"]
    # GMRES settings ----------------------------
    Inner_Iter = GMRES_settings.Inner_Iter
    Outer_Iter = GMRES_settings.Outer_Iter
    # -------------------------------------------
    m = size(incidence_selection["A"], 1)
    n = size(incidence_selection["A"], 2)
    ns = size(incidence_selection["Gamma"], 2)
    w = 2 .* pi .* freq
    nfreq = length(w)
    is = zeros(n, 1)
    S = zeros(ComplexF64, size(ports.port_nodes, 1), size(ports.port_nodes, 1), length(freq))
    Vrest = zeros(ComplexF64, m + n + ns, size(ports.port_nodes, 1))
    invP = sparse(1:ns, 1:ns, 1 ./ diagonals["P"],ns,ns)
    R_chiusura = 50

    # YleArray = []
    # ZSelfArray = []
    # LArray = []
    # UArray = []
    # PArray = []
    # QArray = []
    # SSArray = []
    
    for k = 1:nfreq
        println("Freq n=$k - Freq Tot=$nfreq")
        # Questa parte la vedremo in futuro (ignoratela per ora tanto QS_Rcc_FW è impostato ad 1)-------------------------------------------------------
        if QS_Rcc_FW == 2
            FFTCLp_rebuilted = compute_Circulant_Lp_Rcc(FFTCLp, escalings, freq[k] ./ escalings["freq"])
            FFTCP_rebuilted = compute_Circulant_P_sup_Rcc(FFTCP, escalings, freq[k] ./ escalings["freq"])
        elseif QS_Rcc_FW == 3
            FFTCLp_rebuilted = compute_Circulant_Lp_FW(FFTCLp, escalings, freq[k] ./ escalings["freq"])
            FFTCP_rebuilted = compute_Circulant_P_sup_FW(FFTCP, escalings, freq[k] ./ escalings["freq"])
        end
        # ---------------------------------------------------------------------------

        Yle = build_Yle_S(lumped_elements, [], ports, escalings, n, w[k] / escalings["freq"], R_chiusura)
        #push!(YleArray, Yle)
        Z_self = compute_Z_self(diagonals["R"], diagonals["Cd"], w[k])
        #push!(ZSelfArray, Z_self)
        Zs = escalings["R"] * (Zs_info["Zs"] * sqrt(w[k] / escalings["freq"]))
        ind_to_put_zero_Z_self = findall((real.(Zs[Zs_info["surface_edges"]]) .- real.(Z_self[Zs_info["surface_edges"]])) .> 0)
        ind_to_put_zero_Zs = findall((real.(Zs[Zs_info["surface_edges"]]) .- real.(Z_self[Zs_info["surface_edges"]])) .< 0)
        Z_self[Zs_info["surface_edges"][ind_to_put_zero_Z_self]] .= 0 .+ 1im * imag.(Z_self[Zs_info["surface_edges"][ind_to_put_zero_Z_self]])
        Zs[Zs_info["surface_edges"][ind_to_put_zero_Zs]] .= 0 .+ 1im * imag.(Zs[Zs_info["surface_edges"][ind_to_put_zero_Zs]])
        DZ = Z_self .+ real.(Zs)
        DZ .= DZ .+ 1im * w[k] * diagonals["Lp"]
        invZ = sparse(1:m, 1:m, 1 ./ DZ[:],m,m)
        
        # --------------------- preconditioner ------------------------
        SS = Yle+(prod_real_complex(transpose(incidence_selection["A"]) , prod_complex_real(invZ , incidence_selection["A"])) + 1im * w[k] * incidence_selection["Gamma"] * invP * transpose(incidence_selection["Gamma"]))
        @time F = LinearAlgebra.lu(SS)
        # push!(LArray, F.L)
        # push!(UArray, F.U)
        # push!(PArray, F.p)
        # push!(QArray, F.q)
        # push!(SSArray, SS)
        
    
        # --------------------------------------------------------------
        for c1 = 1:size(ports.port_nodes, 1)
            n1 = convert(Int64,ports.port_nodes[c1, 1])
            n2 = convert(Int64,ports.port_nodes[c1, 2])
            is[n1] = 1 * escalings["Is"]
            is[n2] = -1 * escalings["Is"]
            tn = precond_3_3_Kt(F, invZ, invP, incidence_selection["A"], incidence_selection["Gamma"], m, ns, is)
    
            products_law = x -> ComputeMatrixVector(x, w[k], incidence_selection, FFTCP, FFTCLp, DZ, Yle, expansions, invZ, invP, F)
            prodts = LinearMap{ComplexF64}(products_law, n + m + ns, n + m + ns)
            x0::Vector{ComplexF64}=Vrest[:, c1];
            if QS_Rcc_FW == 1
                #V, flag, relres, iter, resvec = gmres(ComputeMatrixVector, tn, Inner_Iter, GMRES_settings.tol[k], Outer_Iter, [], [], Vrest[:, c1], w[k], incidence_selection, FFTCP, FFTCLp, DZ, Yle, expansions, invZ, invP, L1, U1, P1, Q1)
                V, info = gmres!(x0,prodts, tn , reltol=GMRES_settings.tol[k], restart=Inner_Iter, maxiter=Inner_Iter, initially_zero=false, log=true, verbose=false)
            else
                #V, flag, relres, iter, resvec = gmres(ComputeMatrixVector, tn, Inner_Iter, GMRES_settings.tol[k], Outer_Iter, [], [], Vrest[:, c1], w[k], incidence_selection, FFTCP_rebuilted, FFTCLp_rebuilted, DZ, Yle, expansions, invZ, invP, L1, U1, P1, Q1)
            end
            
            println(info)
            Vrest[:, c1] = V
            is[n1] = 0
            is[n2] = 0
            for c2 = c1:size(ports.port_nodes, 1)
                n3 = convert(Int64, ports.port_nodes[c2, 1])
                n4 = convert(Int64, ports.port_nodes[c2, 2])
                if c1 == c2
                    S[c1, c2, k] = (2 * (V[m+ns+n3] - V[m+ns+n4]) - R_chiusura) / R_chiusura
                else
                    S[c1, c2, k] = (2 * (V[m+ns+n3] - V[m+ns+n4])) / R_chiusura
                end
                S[c2, c1, k] = S[c1, c2, k]
            end
        end
    end
    # matwrite("/Users/edgardovittoria/Downloads/FFTSolverInput.mat", Dict(
    #         "ZSelfArray" => ZSelfArray,
    #         "YleArray" => YleArray,
    #         "LArray" => LArray,
    #         "UArray" => UArray,
    #         "PArray" => PArray,
    #         "QArray" => QArray,
    #         "SSArray" => SSArray
    #     ))
    out = Dict()
    out["S"] = S
    out["Z"] = s2z(S, R_chiusura)
    out["Y"] = s2y(S, R_chiusura)
    out["f"] = freq ./ escalings["freq"]
    return out
end

function ComputeMatrixVector(x_in, w, incidence_selection, FFTCP, FFTCLp, DZ, Yle, expansions, invZ, invP, lu)
    x::Vector{ComplexF64}=x_in[:,1];
    m = size(incidence_selection["A"], 1)
    ns = size(incidence_selection["Gamma"], 2)
    I = x[1:m]
    Q = x[m+1:m+ns]
    Phi = x[m+ns+1:end]
    # Lp * I ---------------------------------------------------------------
    mx = incidence_selection["mx"]
    my = incidence_selection["my"]
    mz = incidence_selection["mz"]
    Y1 = zeros(ComplexF64, m)
    ind_aux_Lp = Dict{Int,Vector{Int}}()
    ind_aux_Lp[1] = 1:mx
    ind_aux_Lp[2] = mx+1:mx+my
    ind_aux_Lp[3] = mx+my+1:mx+my+mz
    for cont = 1:3
        Nx = size(FFTCLp[cont, 1], 1) ÷ 2
        Ny = size(FFTCLp[cont, 1], 2) ÷ 2
        Nz = size(FFTCLp[cont, 1], 3) ÷ 2
        Ired = I[ind_aux_Lp[cont]]
        I_exp = prod_real_complex(expansions["mat_map_Lp"][cont, 1] , Ired)
        CircKT = reshape(I_exp, Nx, Ny, Nz)
        padded_CircKt = zeros(ComplexF64, 2*Nx,2*Ny,2*Nz)
        padded_CircKt[1:size(CircKT,1), 1:size(CircKT,2), 1:size(CircKT,3)] = CircKT
        Chi = ifft(FFTCLp[cont, 1] .* fft(padded_CircKt))
        Y1[ind_aux_Lp[cont]] = Y1[ind_aux_Lp[cont]] + prod_real_complex(transpose(expansions["mat_map_Lp"][cont, 1]) , reshape(Chi[1:Nx, 1:Ny, 1:Nz], Nx * Ny * Nz, 1))
    end
    Y1 = 1im * w * Y1 + DZ .* I + prod_real_complex(incidence_selection["A"] , Phi)
    # ---------------- P * Q ---------------------------------------------
    Y2 = zeros(ns)
    for cont1 = 1:3
        for cont2 = cont1:3
            Nx = size(FFTCP[cont1, cont2], 1) ÷ 2
            Ny = size(FFTCP[cont1, cont2], 2) ÷ 2
            Nz = size(FFTCP[cont1, cont2], 3) ÷ 2
            Q_exp = prod_real_complex(expansions["exp_P"][cont1, cont2] , Q)
            CircKT = reshape(Q_exp, Nx, Ny, Nz)
            padded_CircKt = zeros(ComplexF64, 2*Nx,2*Ny,2*Nz)
            padded_CircKt[1:size(CircKT,1), 1:size(CircKT,2), 1:size(CircKT,3)] = CircKT
            Chi = ifft(FFTCP[cont1, cont2] .* fft(padded_CircKt))
            Y2 = Y2 + transpose(expansions["exp_P"][cont2, cont1]) * (reshape(Chi[1:Nx, 1:Ny, 1:Nz], Nx * Ny * Nz, 1))
            if cont1 != cont2
                Q_exp = prod_real_complex(expansions["exp_P"][cont2, cont1] , Q)
                CircKT = reshape(Q_exp, Nx, Ny, Nz)
                padded_CircKt = zeros(ComplexF64, 2*Nx,2*Ny,2*Nz)
                padded_CircKt[1:size(CircKT,1), 1:size(CircKT,2), 1:size(CircKT,3)] = CircKT
                Chi = ifft(FFTCP[cont1, cont2] .* fft(padded_CircKt))
                Y2 = Y2 + prod_real_complex(transpose(expansions["exp_P"][cont1, cont2]) , (reshape(Chi[1:Nx, 1:Ny, 1:Nz], Nx * Ny * Nz, 1)))
            end
        end
    end
    Y2 = Y2 - prod_real_complex(transpose(incidence_selection["Gamma"]) , Phi)
    Y3 = -1.0*(prod_real_complex(transpose(incidence_selection["A"]) , I)) + prod_real_complex(Yle , Phi) + 1im * w * (prod_real_complex(incidence_selection["Gamma"] , Q))
    MatrixVector = precond_3_3_vector(lu, invZ, invP, incidence_selection["A"], incidence_selection["Gamma"], w, Y1, Y2, Y3)
    return MatrixVector
end

# function precond_3_3_Kt(L1, U1, P1, Q1, invZ, invP, incidence_selection, n1, n2, X3)
#     n3 = length(X3)
#     i1 = 1:n1
#     i2 = n1+1:n1+n2
#     i3 = n1+n2+1:n1+n2+n3
#     Y = zeros(n1 + n2 + n3, 1)
#     M5 = Q1 * (U1 \ (L1 \ (P1 * X3)))
#     Y[i1] = Y[i1] - invZ * (incidence_selection["A"] * M5)
#     Y[i2] = Y[i2] + invP * (transpose(incidence_selection["Gamma"]) * M5)
#     Y[i3] = Y[i3] + M5
#     return Y
# end

function precond_3_3_Kt(F, invZ, invP, A,Gamma, n1,n2, X3)

    n3 = length(X3)

    i1 = range(1, stop=n1)
    i2 = range(n1+1, stop=n1 + n2)
    i3 = range(n1 + n2 + 1, stop=n1 + n2 + n3)

    Y = zeros(ComplexF64, n1 + n2 + n3)

    #println(X3)

    M5 = F\X3
    # y = F.L \ (F.Rs .* X3)[F.p]
    # z = F.U \ y
    # M5 = z[invperm(F.q)]
    #display(M5)

    Y[i1] .= Y[i1] .- 1.0*(prod_real_complex(invZ, prod_real_complex(A, M5)))
    Y[i2] .= Y[i2] .+ (prod_real_complex(invP, prod_real_complex(transpose(Gamma), M5)))
    Y[i3] .= Y[i3] .+ M5

    return Y
end

# function precond_3_3_vector(L1,U1,P1,Q1,invZ,invP,incidence_selection,w,X1,X2,X3)
#     n1=length(X1)
#     n2=length(X2)
#     n3=length(X3)
#     i1=1:n1
#     i2=n1+1:n1+n2
#     i3=n1+n2+1:n1+n2+n3
#     Y=zeros(n1+n2+n3,1)
#     M1 = invZ*X1
#     M2 = (Q1*(U1\(L1\(P1*(transpose(incidence_selection["A"])*M1)))))
#     M3 = invP*X2
#     M4 = (Q1*(U1\(L1\(P1*(incidence_selection["Gamma"]*M3)))))
#     M5 = Q1*(U1\(L1\(P1*X3)))
#     Y[i1] .= Y[i1] .+ invZ*X1 .- invZ*(incidence_selection["A"]*M2)
#     Y[i1] .= Y[i1] .+ 1im*w*(invZ*(incidence_selection["A"]*M4))
#     Y[i1] .= Y[i1] .- invZ*(incidence_selection["A"]*M5)
#     Y[i2] .= Y[i2] .+ invP*(transpose(incidence_selection["Gamma"])*M2)
#     Y[i2] .= Y[i2] .+ invP*X2 .- 1im*w*invP*(transpose(incidence_selection["Gamma"])*M4)
#     Y[i2] .= Y[i2] .+ invP*(transpose(incidence_selection["Gamma"])*M5)
#     Y[i3] .= Y[i3] .+ M2
#     Y[i3] .= Y[i3] .- 1im*w*M4
#     Y[i3] .= Y[i3] .+ M5
#     return Y
# end
function precond_3_3_vector(F,invZ,invP,A,Gamma,w,X1,X2,X3)

    n1=length(X1)
    n2=length(X2)
    n3=length(X3)

    i1=range(1, stop=n1)
    i2=range(n1+1,stop=n1+n2)
    i3=range(n1+n2+1,stop=n1+n2+n3)

    Y=zeros(ComplexF64 , n1+n2+n3)


    #da rivedere
    M1 = prod_real_complex(invZ, X1)
    #M1 = csc_matrix.*(invZ, X1)
    M2 = F\(prod_real_complex(transpose(A), M1))
    # y = F.L \ (F.Rs .* prod_real_complex(transpose(A), M1))[F.p]
    # z = F.U \ y
    # M2 = z[invperm(F.q)]
    #M2 = lu_S.solve(csc_matrix.dot(A.transpose(),M1))
    M3 = prod_real_complex(invP, X2)
    #M3 = csc_matrix.*(invP,X2)
    M4 = F\prod_real_complex(Gamma, M3)
    # y = F.L \ (F.Rs .* prod_real_complex(Gamma, M3))[F.p]
    # z = F.U \ y
    # M4 = z[invperm(F.q)]
    #M4 = lu_S.solve(csc_matrix.*(Gamma,M3))
    M5 = F\X3
    # y = F.L \ (F.Rs .* X3)[F.p]
    # z = F.U \ y
    # M5 = z[invperm(F.q)]

    # for i in i1
    #     Y[i] = Y[i]+M1-1.0*(*((invZ),*((A), M2)))
    #     Y[i] = Y[i]+1im*w*(*((invZ),*((A), M4)))
    #     Y[i] = Y[i]-1.0*(*((invZ),*((A), M5)))
    # end
    Y[i1] .= Y[i1] .+ M1-1.0*(prod_real_complex((invZ),prod_real_complex((A), M2)))
    Y[i1] .= Y[i1] .+ 1im*w*(prod_real_complex((invZ),prod_real_complex((A), M4)))
    Y[i1] .= Y[i1] .- 1.0*(prod_real_complex((invZ),prod_real_complex((A), M5)))

    #Y[np.ix_(i1)] = Y[np.ix_(i1)]+M1-1.0*csc_matrix.*(invZ,csc_matrix.*(A,M2))
    #Y[np.ix_(i1)] = Y[np.ix_(i1)]+1im*w*csc_matrix.*(invZ, csc_matrix.*(A,M4))
    #Y[np.ix_(i1)] = Y[np.ix_(i1)]-1.0*csc_matrix.*(invZ, csc_matrix.*(A, M5))

    
    Y[i2] .= Y[i2] .+ (prod_real_complex(invP,prod_real_complex((transpose(Gamma)), M2)))
    Y[i2] .= Y[i2] .+ M3 - 1im*w*(prod_real_complex(invP,prod_real_complex(transpose(Gamma), M4)))
    Y[i2] .= Y[i2] .+ (prod_real_complex(invP,prod_real_complex(transpose(Gamma), M5)))
    

    # Y[np.ix_(i2)] = Y[np.ix_(i2)]+csc_matrix.*(invP, csc_matrix.*(Gamma.transpose(), M2))
    # Y[np.ix_(i2)] = Y[np.ix_(i2)] + M3 -1im*w*csc_matrix.*(invP, csc_matrix.*(Gamma.transpose(), M4))
    # Y[np.ix_(i2)] = Y[np.ix_(i2)]+csc_matrix.*(invP, csc_matrix.*(Gamma.transpose(), M5))

    
    Y[i3] .= Y[i3] .+ M2
    Y[i3] .= Y[i3] .- 1im*w*M4
    Y[i3] .= Y[i3] .+ M5
    

    # Y[np.ix_(i3)] = Y[np.ix_(i3)]+M2
    # Y[np.ix_(i3)] = Y[np.ix_(i3)]-1im*w*M4
    # Y[np.ix_(i3)] = Y[np.ix_(i3)]+M5

    Y=convert(Array{Complex{Float64}}, Y)
    return Y
end

function prod_real_complex(A,x)
    # A is a N x N real matrix and x is a complex matrix

    N=size(A,1);
    y=zeros(ComplexF64 , N, 1)
    y=*(A,real.(x))+1im * *(A,imag.(x))
    return y
end

function prod_complex_real(A,x)
    # A is a N x N complex matrix and x is a real matrix

    N=size(A,1);
    y=zeros(ComplexF64 , N, 1)
    y=*(real.(A),x)+1im * *(imag.(A),x)
    return y
end

function s2z(S,Zo)
    num_ports=size(S)[1]
    nfreq=size(S)[3]
    Z = zeros(ComplexF64 , num_ports, num_ports, nfreq)
    Id = Matrix{Int64}(I, num_ports, num_ports)
    for cont in range(1, stop=nfreq)
        Z[:,:,cont]=Zo*((Id-1.0*S[:,:,cont])\(Id+S[:,:,cont]))
    
    end
    return Z
end

function s2y(S,Zo)
    num_ports=size(S)[1]
    nfreq=size(S)[3]
    Y = zeros(ComplexF64 , num_ports, num_ports, nfreq)
    Id = Matrix{Int64}(I, num_ports, num_ports)
    for cont in range(1, stop=nfreq)
        Y[:,:,cont]=Zo*((Id+S[:,:,cont])\(Id-1.0*S[:,:,cont]))
    end
    return Y
end
