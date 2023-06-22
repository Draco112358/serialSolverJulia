include("From_3D_to_1D.jl")
using Base.Sort

function create_expansion_ind_Lp_x_grids_v2(grids::Array{Array{Int64,3},1}, map::Array{Int64,2}, l_ind::Int64, mapping_Vox::Array{Int64,1}, nodes::Array{Float64,1}, nodes_red::Array{Float64,1})
    num_grids = length(grids)
    Nx, Ny, Nz = size(grids[1])
    Ind_out = zeros(Int64, l_ind, 2)
    pos = 0
    for cont3 in 1:Nz
        for cont2 in 1:Ny
            for cont in 1:Nx-1
                for k in 1:num_grids
                    if grids[k][cont,cont2,cont3] && grids[k][cont+1,cont2,cont3]
                        nn1 = searchsortedfirst(nodes_red, nodes[mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]])
                        nn2 = searchsortedfirst(nodes_red, nodes[mapping_Vox[From_3D_to_1D(cont+1, cont2, cont3, Nx, Ny)]])
                        if abs(nn1-nn2) > 1e-8
                            pos += 1
                            Ind_out[pos,1] = From_3D_to_1D(cont, cont2, cont3, Nx-1, Ny)
                            Ind_out[pos,2] = map[From_3D_to_1D(cont, cont2, cont3, Nx-1, Ny)]
                            break
                        end
                    end
                end
            end
        end
    end
    return Ind_out
end

function bin_search(num::Float64, A::Array{Float64,1})
    index = 0
    n = length(A)
    left = 1
    right = n
    while left <= right
        mid = ceil(Int64, (left + right) / 2)
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
