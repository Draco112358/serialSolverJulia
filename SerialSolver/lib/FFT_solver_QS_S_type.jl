include("build_Yle_S.jl")
include("compute_Z_self.jl")
include("gmres_custom.jl")
using SparseArrays, SuperLU, IterativeSolvers, FFTW, FourierTools, MKL
using LinearMaps, MAT, LinearAlgebra, Krylov, KrylovMethods, ILUZero,FLoops


function FFT_solver_QS_S_type(freq, escalings, incidence_selection, FFTCP, FFTCLp, diagonals, ports, lumped_elements, expansions, GMRES_settings, Zs_info, QS_Rcc_FW)
    FFTW.set_num_threads(12)
    BLAS.set_num_threads(1)
    # println(BLAS.get_num_threads())
    #
    #println("qui")
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
    is = zeros(Float64, n, 1)
    S = zeros(ComplexF64, size(ports.port_nodes, 1), size(ports.port_nodes, 1), length(freq))
    Vrest = zeros(ComplexF64, m + n + ns, size(ports.port_nodes, 1))
    invP = sparse(1:ns, 1:ns, 1 ./ diagonals["P"],ns,ns)
    R_chiusura = 50.0
    PVector:: Vector{Any} = []
    PLIVector:: Vector{Any} = []
    for cont = 1:3
            Nx = size(FFTCLp[cont, 1], 1) ÷ 2
            Ny = size(FFTCLp[cont, 1], 2) ÷ 2
            Nz = size(FFTCLp[cont, 1], 3) ÷ 2
            padded_CircKt = zeros(ComplexF64, 2*Nx,2*Ny,2*Nz)
            P = plan_fft(padded_CircKt, flags=FFTW.MEASURE)
            push!(PVector, P)
            push!(PLIVector, plan_ifft(FFTCLp[cont, 1] .* (P*padded_CircKt), flags=FFTW.MEASURE))
    end
    P2Vector:: Matrix{Any} = zeros(3, 3)
    PLI2Vector:: Matrix{Any} = zeros(3, 3)
    for cont1 = 1:3
        for cont2 = cont1:3
            Nx = size(FFTCP[cont1, cont2], 1) ÷ 2
            Ny = size(FFTCP[cont1, cont2], 2) ÷ 2
            Nz = size(FFTCP[cont1, cont2], 3) ÷ 2
            padded_CircKt = zeros(ComplexF64, 2*Nx,2*Ny,2*Nz)
            #Chi = ifft(FFTCP[cont1, cont2] .* fft(padded_CircKt))
            P = plan_fft(padded_CircKt, flags=FFTW.MEASURE)
            P2Vector[cont1, cont2] = P
            PLI2Vector[cont1, cont2] = plan_ifft(FFTCP[cont1, cont2] .* (P*padded_CircKt), flags=FFTW.MEASURE)
        end
    end
    
    for k = 1:nfreq
        #println("Freq n=$k - Freq Tot=$nfreq")
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
        Z_self = compute_Z_self(diagonals["R"], diagonals["Cd"], w[k])
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
        #SS = mxarray(SS)
        # L1 = @mget L
        # U1 = @mget U
        # P1 = @mget P
        # Q1 = @mget Q

        #F = @time splu(SS)
        
        
        F = lu(SS)
        #F = SS
        #F = @time ilu0(SS)
        #ps = PardisoSolver()
        #set_iparm!(ps,1,1)
        #set_iparm!(ps, 8, 10)
        #set_matrixtype!(ps, Pardiso.COMPLEX_NONSYM)
        #set_solver!(ps, Pardiso.ITERATIVE_SOLVER)
        
        # --------------------------------------------------------------
        for c1 = 1:size(ports.port_nodes, 1)
            n1 = convert(Int64,ports.port_nodes[c1, 1])
            n2 = convert(Int64,ports.port_nodes[c1, 2])
            is[n1] = 1 * escalings["Is"]
            is[n2] = -1 * escalings["Is"]
            tn = precond_3_3_Kt(F, invZ, invP, incidence_selection["A"], incidence_selection["Gamma"], m, ns, is)
    
            # products_law = x ->   ComputeMatrixVector(x, w[k], incidence_selection, FFTCP, FFTCLp, DZ, Yle, expansions, invZ, invP, F)
            # prodts = LinearMap{ComplexF64}(products_law, n + m + ns, n + m + ns)
            # x0::Vector{ComplexF64}=Vrest[:, c1];
        
            if QS_Rcc_FW == 1
                V, flag, relres, iter, resvec = gmres_custom(tn, false, GMRES_settings.tol[k], Inner_Iter, Vrest[:, c1], w[k], incidence_selection, FFTCP, FFTCLp, DZ, Yle, expansions, invZ, invP, F, PLIVector, PVector, PLI2Vector, P2Vector)
                tot_iter_number = (iter[1] - 1) * Inner_Iter + iter[2] + 1
                if (flag == 0)
                    println("Flag $flag - Iteration = $k - Convergence reached, number of iterations:$tot_iter_number")
                end

                if (flag == 1)
                    println("Flag $flag - Iteration = $k - Convergence not reached, number of iterations:$Inner_Iter")
                end
                #V, info = IterativeSolvers.gmres!(x0, prodts, tn; reltol=GMRES_settings.tol[k], restart=Inner_Iter, maxiter=Inner_Iter, initially_zero=false, log=true, verbose=true)
                #V, info = IterativeSolvers.gmres(prodts, tn; reltol=GMRES_settings.tol[k], restart=Inner_Iter, maxiter=Inner_Iter, initially_zero=false, log=true, verbose=true)
                #println(info)
                
                # (V, stats) = Krylov.gmres(prodts,tn,x0;restart=false, memory=Outer_Iter, reorthogonalization=false,rtol=GMRES_settings.tol[k], itmax=Inner_Iter,verbose=1, history=false)
            
                # if (stats.solved==true)        
                #     println("convergence reached, number of iterations: "*string(stats.niter))
                # else
                #     println("convergence not reached, number of iterations: "*string(stats.niter))
                # end

                #V,flag,err,iter,resvec = KrylovMethods.gmres(prodts,tn,Inner_Iter,tol=GMRES_settings.tol[k],maxIter=1,x=x0,out=1)

                # if (flag==0)        
                #     println("convergence reached, number of iterations: "*string(length(resvec)))
                # elseif(flag==-1)
                #     println("convergence not reached, number of iterations: "*string(length(resvec)))
                # else
                #     println("RHS all zeros, number of iterations: "*string(length(resvec)))
                # end
            else
                #V, flag, relres, iter, resvec = gmres(ComputeMatrixVector, tn, Inner_Iter, GMRES_settings.tol[k], Outer_Iter, [], [], Vrest[:, c1], w[k], incidence_selection, FFTCP_rebuilted, FFTCLp_rebuilted, DZ, Yle, expansions, invZ, invP, L1, U1, P1, Q1)
            end
            
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
    out = Dict()
    out["S"] = S
    out["Z"] = s2z(S, R_chiusura)
    out["Y"] = s2y(S, R_chiusura)
    out["f"] = freq ./ escalings["freq"]
    return out
end

function ComputeMatrixVector(x, w, incidence_selection, FFTCP, FFTCLp, DZ, Yle, expansions, invZ, invP, lu, PLIVector, PVector, PLI2Vector, P2Vector)
    # @profile begin
        
    # end
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
            #mat"Chi=ifftn($FFTCLp{$cont,1}.*fftn($CircKT,[2*$Nx,2*$Ny,2*$Nz]))"
            # P = plan_fft(padded_CircKt, flags=FFTW.MEASURE)
            # PLI = plan_ifft(FFTCLp[cont, 1] .* (P*padded_CircKt), flags=FFTW.MEASURE)
            # Chi = PLI * (FFTCLp[cont, 1] .* (P*padded_CircKt))
            Chi = customIfftOptimized(PLIVector, PVector, padded_CircKt, FFTCLp, cont)
            #Chi = customIfft(padded_CircKt, FFTCLp, cont)
            #Chi = ifft(FFTCLp[cont, 1] .* fft(padded_CircKt))
            #Chi = @mget Chi
            Y1[ind_aux_Lp[cont]] = Y1[ind_aux_Lp[cont]] + prod_real_complex(transpose(expansions["mat_map_Lp"][cont, 1]) , reshape(Chi[1:Nx, 1:Ny, 1:Nz], Nx * Ny * Nz, 1))
    end
    Y1 = 1im * w * Y1 + DZ .* I + prod_real_complex(incidence_selection["A"] , Phi)
    # ---------------- P * Q ---------------------------------------------
    Y2 = zeros(ComplexF64,ns)
    for cont1 = 1:3
        for cont2 = cont1:3
            Nx = size(FFTCP[cont1, cont2], 1) ÷ 2
            Ny = size(FFTCP[cont1, cont2], 2) ÷ 2
            Nz = size(FFTCP[cont1, cont2], 3) ÷ 2
            Q_exp = prod_real_complex(expansions["exp_P"][cont1, cont2] , Q)
            CircKT = reshape(Q_exp, Nx, Ny, Nz)
            padded_CircKt = zeros(ComplexF64, 2*Nx,2*Ny,2*Nz)
            padded_CircKt[1:size(CircKT,1), 1:size(CircKT,2), 1:size(CircKT,3)] = CircKT
            #Chi = ifft(FFTCP[cont1, cont2] .* fft(padded_CircKt))
            # P = plan_fft(padded_CircKt, flags=FFTW.MEASURE)
            # PLI = plan_ifft(FFTCP[cont1, cont2] .* (P*padded_CircKt), flags=FFTW.MEASURE)
            # Chi = PLI * (FFTCP[cont1, cont2] .* (P*padded_CircKt))
            #Chi = @time customIfft2(padded_CircKt, FFTCP, cont1, cont2)
            Chi = customIfftOptimized2(PLI2Vector, P2Vector, padded_CircKt, FFTCP, cont1, cont2)
            
            # mat"Chi=ifftn($FFTCP{$cont1,$cont2}.*fftn($CircKT,[2*$Nx,2*$Ny,2*$Nz]))"
            # Chi = @mget Chi
            Y2 = Y2 + prod_real_complex(transpose(expansions["exp_P"][cont2, cont1]) , (reshape(Chi[1:Nx, 1:Ny, 1:Nz], Nx * Ny * Nz, 1)))
            if cont1 != cont2
                Q_exp = prod_real_complex(expansions["exp_P"][cont2, cont1] , Q)
                CircKT = reshape(Q_exp, Nx, Ny, Nz)
                padded_CircKt = zeros(ComplexF64, 2*Nx,2*Ny,2*Nz)
                padded_CircKt[1:size(CircKT,1), 1:size(CircKT,2), 1:size(CircKT,3)] = CircKT
                #Chi = ifft(FFTCP[cont1, cont2] .* fft(padded_CircKt))
                # P = plan_fft(padded_CircKt)
                # PLI = plan_ifft(FFTCP[cont1, cont2] .* (P*padded_CircKt))
                #Chi = PLI * (FFTCP[cont1, cont2] .* (P*padded_CircKt))
                #Chi = @time customIfft2(padded_CircKt, FFTCP, cont1, cont2)
                Chi = customIfftOptimized2(PLI2Vector, P2Vector, padded_CircKt, FFTCP, cont1, cont2)
                # mat"Chi=ifftn($FFTCP{$cont1,$cont2}.*fftn($CircKT,[2*$Nx,2*$Ny,2*$Nz]))"
                # Chi = @mget Chi
                Y2 = Y2 + prod_real_complex(transpose(expansions["exp_P"][cont1, cont2]) , (reshape(Chi[1:Nx, 1:Ny, 1:Nz], Nx * Ny * Nz, 1)))
            end
        end
    end
    Y2 = Y2 - prod_real_complex(transpose(incidence_selection["Gamma"]) , Phi)
    #Y2 = Y2 - threaded_mul!([], transpose(incidence_selection["Gamma"]), Phi)
    Y3 = -1.0*(prod_real_complex(transpose(incidence_selection["A"]) , I)) + prod_real_complex(Yle , Phi) + 1im * w * (prod_real_complex(incidence_selection["Gamma"] , Q))
    #Y3 = -1.0*(prod_real_complex(transpose(incidence_selection["A"]) , I)) + threaded_mul!([], Yle, Phi) + 1im * w * (prod_real_complex(incidence_selection["Gamma"] , Q))
    MatrixVector = precond_3_3_vector(lu, invZ, invP, incidence_selection["A"], incidence_selection["Gamma"], w, Y1, Y2, Y3)
    return MatrixVector
end


function precond_3_3_Kt(F, invZ, invP, A,Gamma, n1,n2, X3)
    n3 = length(X3)

    i1 = range(1, stop=n1)
    i2 = range(n1+1, stop=n1 + n2)
    i3 = range(n1 + n2 + 1, stop=n1 + n2 + n3)

    Y = zeros(ComplexF64, n1 + n2 + n3, 1)

    #M5 = Q1*(U1\(L1\(P1*X3)));
    M5 = F\X3
    #M5 = zeros(ComplexF64, size(X3)[1], size(X3)[2])
    
    #Pardiso.solve!(ps, M5, F, convert(Matrix{ComplexF64}, X3))

    Y[i1] .= Y[i1] .- 1.0*(prod_real_complex(invZ, prod_real_complex(A, M5)))
    Y[i2] .= Y[i2] .+ (prod_real_complex(invP, prod_real_complex(transpose(Gamma), M5)))
    Y[i3] .= Y[i3] .+ M5


    return Y
    
end

function precond_3_3_vector(F,invZ,invP,A,Gamma,w,X1,X2,X3)
    
    n1=length(X1)
    n2=length(X2)
    n3=length(X3)

    i1=range(1, stop=n1)
    i2=range(n1+1,stop=n1+n2)
    i3=range(n1+n2+1,stop=n1+n2+n3)

    Y=zeros(ComplexF64 , n1+n2+n3)
    M1 = prod_real_complex(invZ, X1)
   # M2 = (Q1*(U1\(L1\(P1*(prod_real_complex(transpose(A), M1))))))
    M2 = F\(prod_real_complex(transpose(A), M1))  
    #M2 = solve(ps, F, prod_real_complex(transpose(A), M1))
    M3 = prod_real_complex(invP, X2) 
    #M4 = (Q1*(U1\(L1\(P1*(prod_real_complex(Gamma, M3))))));
    #M5 = Q1*(U1\(L1\(P1*X3)));  
    M4 = F\prod_real_complex(Gamma, M3) 
    #M4 = solve(ps, F, prod_real_complex(Gamma, M3))
    M5 = F\X3
    #M5 = solve(ps, F, X3)
    
    Y[i1] .= Y[i1] .+ M1-1.0*(prod_real_complex((invZ),prod_real_complex((A), M2)))
    Y[i1] .= Y[i1] .+ 1im*w*(prod_real_complex((invZ),prod_real_complex((A), M4)))
    Y[i1] .= Y[i1] .- 1.0*(prod_real_complex((invZ),prod_real_complex((A), M5)))
    
    Y[i2] .= Y[i2] .+ (prod_real_complex(invP,prod_real_complex((transpose(Gamma)), M2)))
    Y[i2] .= Y[i2] .+ M3 - 1im*w*(prod_real_complex(invP,prod_real_complex(transpose(Gamma), M4)))
    Y[i2] .= Y[i2] .+ (prod_real_complex(invP,prod_real_complex(transpose(Gamma), M5)))
    
    Y[i3] .= Y[i3] .+ M2
    Y[i3] .= Y[i3] .- 1im*w*M4
    Y[i3] .= Y[i3] .+ M5
    
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

function customIfftOptimized(PLIVector, PVector, padded_CircKt, FFTCLp, cont)
    return PLIVector[cont] * (FFTCLp[cont, 1] .* (PVector[cont]*padded_CircKt))
end

function customIfftOptimized2(PLI2Vector, P2Vector, padded_CircKt, FFTCP, cont1, cont2)
    return PLI2Vector[cont1,cont2] * (FFTCP[cont1, cont2] .* (P2Vector[cont1,cont2]*padded_CircKt))
end


function customIfft(padded_CircKt, FFTCLp, cont)
    # P = plan_fft(padded_CircKt, flags=FFTW.MEASURE)
    # y = deepcopy(padded_CircKt)
    # mul!(y, P, padded_CircKt)
    # PLI = plan_ifft(FFTCLp[cont, 1] .* y, flags=FFTW.MEASURE)
    #PLI = plan_ifft(FFTCLp[cont, 1] .* (P*padded_CircKt), flags=FFTW.MEASURE)
    # Chi = deepcopy(FFTCLp[cont, 1] .* y)
    # mul!(Chi, PLI, (FFTCLp[cont, 1] .* y))
    #Chi = PLI * (FFTCLp[cont, 1] .* (P*padded_CircKt))
    Chi = ifft(FFTCLp[cont, 1] .* fft(padded_CircKt))
    return Chi
end

function customIfft2(padded_CircKt, FFTCP, cont1, cont2)
    # P = plan_fft(padded_CircKt)
    # PLI = plan_ifft(FFTCP[cont1, cont2] .* (P*padded_CircKt))
    # Chi = PLI * (FFTCP[cont1, cont2] .* (P*padded_CircKt))
    Chi = ifft(FFTCP[cont1, cont2] .* fft(padded_CircKt))
    return Chi
end
