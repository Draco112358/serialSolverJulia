function centri_vox=build_center_air_Voxels(grids,sx,sy,sz,min_v)

Nx=size(grids{1},1);
Ny=size(grids{1},2);
Nz=size(grids{1},3);

centri_vox=zeros(Nx*Ny*Nz,3);

num_ele=0;
for cont=1:Nx
    for cont2=1:Ny
        for cont3=1:Nz
            
            num_ele=num_ele+1;
            cx=min_v(1)+sx*(cont-1)+sx/2;
            cy=min_v(2)+sy*(cont2-1)+sy/2;
            cz=min_v(3)+sz*(cont3-1)+sz/2;
            centri_vox((From_3D_to_1D(cont, cont2, cont3, Nx, Ny)),:)=[cx cy cz];
           
        end
    end
end


end