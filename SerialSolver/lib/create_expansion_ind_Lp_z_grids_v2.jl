using Base.Sort: searchsortedfirst

function create_expansion_ind_Lp_z_grids_v2(grids, map, l_ind, mapping_Vox, nodes, nodes_red)
    num_grids = length(grids)
    Nx, Ny, Nz = size(grids[1])
    Ind_out = zeros(l_ind, 2)
    pos = 0
    for cont = 1:Nx, cont2 = 1:Ny, cont3 = 1:Nz-1
        for k = 1:num_grids
            if grids[k][cont, cont2, cont3] && grids[k][cont, cont2, cont3+1]
                nn1 = searchsortedfirst(nodes_red, nodes[mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]])
                nn2 = searchsortedfirst(nodes_red, nodes[mapping_Vox[From_3D_to_1D(cont, cont2, cont3+1, Nx, Ny)]])
                if abs(nn1-nn2) > 1e-8
                    pos += 1
                    Ind_out[pos, 1] = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    Ind_out[pos, 2] = map[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]
                    break
                end
            end
        end
    end
    return Ind_out
end

function bin_search(num, A)
    index = searchsortedfirst(A, num)
    if index > length(A) || A[index] != num
        index = 0
    end
    return index
end
