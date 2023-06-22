include("From_3D_to_1D.jl")
using Base.Iterators

function create_mapping_Ax_v2(grids, mapping_Vox, nodes, nodes_red)
    num_grids = length(grids)
    Nx, Ny, Nz = size(grids[1])
    N_max = (Nx-1)*Ny*Nz
    mapping = zeros(N_max)
    num_ele = 0
    for (cont2, cont3, cont) in product(1:Ny, 1:Nz, 1:Nx-1)
        for k in 1:num_grids
            if grids[k][cont, cont2, cont3] && grids[k][cont+1, cont2, cont3]
                nn1 = searchsortedfirst(nodes_red, nodes[mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]])
                nn2 = searchsortedfirst(nodes_red, nodes[mapping_Vox[From_3D_to_1D(cont+1, cont2, cont3, Nx, Ny)]])
                if abs(nn1-nn2) > 1e-8
                    kkey = From_3D_to_1D(cont, cont2, cont3, Nx-1, Ny)
                    if mapping[kkey] == 0
                        num_ele += 1
                        mapping[kkey] = num_ele
                    end
                    break
                end
            end
        end
    end
    return mapping, num_ele
end

function bin_search(num, A)
    index = searchsortedfirst(A, num)
    if index > length(A) || A[index] != num
        index = 0
    end
    return index
end
