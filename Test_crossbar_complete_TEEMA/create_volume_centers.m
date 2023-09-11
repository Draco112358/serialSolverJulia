function [centri_vox,id_mat]=create_volume_centers(grids,map,num_centri,sx,sy,sz)

global min_v;

centri_vox=zeros(num_centri,3);
id_mat=zeros(num_centri,1);

Nx=size(grids{1},1);
Ny=size(grids{1},2);
Nz=size(grids{1},3);

num_grids=length(grids);

for cont=1:Nx
    for cont2=1:Ny
        for cont3=1:Nz
            for k=1:num_grids
                if(grids{k}(cont,cont2,cont3))
                    pos=map(From_3D_to_1D(cont,cont2,cont3, Nx, Ny));
                    cx=0+sx*(cont-1)+sx/2;
                    cy=0+sy*(cont2-1)+sy/2;
                    cz=0+sz*(cont3-1)+sz/2;
                    centri_vox(pos,:)=[cx cy cz];
                    id_mat(pos)=k;
                    break;
                end
            end
        end
    end
end


end

function pos=From_3D_to_1D(i, j, k, M, N) 
	pos = ((k-1) * M * N) + ((j-1) * M) + i;
end