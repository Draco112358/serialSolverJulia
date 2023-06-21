function [mapping,num_ele]=create_volumes_mapping_v2(grids)

num_grids=length(grids);

Nx=size(grids{1},1);
Ny=size(grids{1},2);
Nz=size(grids{1},3);

mapping=zeros(Nx*Ny*Nz,1);

num_ele=0;
for cont=1:Nx
    for cont2=1:Ny
        for cont3=1:Nz
            for k=1:num_grids
                if(grids{k}(cont,cont2,cont3))
                    num_ele=num_ele+1;
                    mapping(From_3D_to_1D(cont, cont2, cont3, Nx, Ny))=num_ele;
                    break;
                end
            end
        end
    end
end

end