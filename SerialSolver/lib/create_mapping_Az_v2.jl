using LinearAlgebra

function create_mapping_Az_v2(grids, mapping_Vox, nodes, nodes_red)
    num_grids = length(grids)
    Nx, Ny, Nz = size(grids[1])
    N_max = Nx * Ny * (Nz - 1)
    mapping = zeros(Int64, N_max)
    num_ele = 0
    for cont = 1:Nx
        for cont2 = 1:Ny
            for cont3 = 1:Nz-1
                for k = 1:num_grids
                    if grids[k][cont, cont2, cont3] && grids[k][cont, cont2, cont3+1]
                        nn1 = searchsorted(nodes_red, nodes[mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]])
                        nn2 = searchsorted(nodes_red, nodes[mapping_Vox[From_3D_to_1D(cont, cont2, cont3+1, Nx, Ny)]])
                        if abs(nn1-nn2) > 1e-8
                            kkey = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                            if mapping[kkey] == 0
                                num_ele += 1
                                mapping[kkey] = num_ele
                            end
                            break
                        end
                    end
                end
            end
        end
    end
    return mapping, num_ele
end

function bin_search(num, A)
    index = 0
    n = length(A)
    left = 1
    right = n
    while left <= right
        mid = ceil((left + right) / 2)
        if A[mid] == num
            index = mid
            break
        else
            if A[mid] > num
                right = mid - 1
            else
                left = mid + 1
            end
        end
    end
    return index
end
