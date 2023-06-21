function Ind_out=create_expansion_ind_Lp_z_grids_v2(grids,map,l_ind,mapping_Vox,nodes,nodes_red)

num_grids=length(grids);

Nx=size(grids{1},1);
Ny=size(grids{1},2);
Nz=size(grids{1},3);

Ind_out=zeros(l_ind,2);

pos=0;
for cont=1:Nx
    for cont2=1:Ny
        for cont3=1:Nz-1
            
            for k=1:num_grids
                if(grids{k}(cont,cont2,cont3)) && (grids{k}(cont,cont2,cont3+1))
                    nn1=bin_search(nodes(mapping_Vox(From_3D_to_1D(cont, cont2, cont3, Nx, Ny))),nodes_red);
                    nn2=bin_search(nodes(mapping_Vox(From_3D_to_1D(cont, cont2, cont3+1, Nx, Ny))),nodes_red);
                    
                    if abs(nn1-nn2)>1e-8
                        pos=pos+1;
                        Ind_out(pos,1)=From_3D_to_1D(cont, cont2, cont3, Nx, Ny);
                        Ind_out(pos,2)=map(From_3D_to_1D(cont, cont2, cont3, Nx, Ny));
                        break;
                    end
                end
            end
        end
    end
end

end


function index = bin_search(num,A)

index=0;
n=length(A);
left = 1;
right = n;

while left <= right
    mid = ceil((left + right) / 2);
    
    if A(mid) == num
        index = mid;
        break;
    else
        if A(mid) > num
            right = mid - 1;
        else
            left = mid + 1;
        end
    end
end

end