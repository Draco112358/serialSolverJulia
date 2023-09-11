function centri_vox=build_center_P_Voxels(min_v,Nx,Ny,Nz,sx,sy,sz)

centri_vox.p12_se=zeros((Nx)*(Ny+1)*(Nz),3);
centri_vox.p34_se=zeros((Nx+1)*(Ny)*(Nz),3);
centri_vox.p56_se=zeros((Nx)*(Ny)*(Nz+1),3);

centri_vox.p1234=zeros((2*(Nx+1)-1)*(2*(Ny+1)-1)*(Nz),3);
centri_vox.p1256=zeros((Nx)*(2*(Ny+1)-1)*(2*(Nz+1)-1),3);
centri_vox.p3456=zeros((2*(Nx+1)-1)*(Ny)*(2*(Nz+1)-1),3);

for cont=1:Nx
    for cont2=1:Ny+1
        for cont3=1:Nz
            
            cx=0+sx*(cont-1)+sx/2;
            cy=0+sy*(cont2-1)+sy/2;
            cz=0+sz*(cont3-1)+sz/2;
            centri_vox.p12_se((From_3D_to_1D(cont, cont2, cont3, Nx, Ny+1)),:)=[cx cy cz];
           
        end
    end
end

for cont=1:Nx+1
    for cont2=1:Ny
        for cont3=1:Nz
            
            cx=0+sx*(cont-1)+sx/2;
            cy=0+sy*(cont2-1)+sy/2;
            cz=0+sz*(cont3-1)+sz/2;
            centri_vox.p34_se((From_3D_to_1D(cont, cont2, cont3, Nx+1, Ny)),:)=[cx cy cz];
           
        end
    end
end

for cont=1:Nx
    for cont2=1:Ny
        for cont3=1:Nz+1
            
            cx=0+sx*(cont-1)+sx/2;
            cy=0+sy*(cont2-1)+sy/2;
            cz=0+sz*(cont3-1)+sz/2;
            centri_vox.p56_se((From_3D_to_1D(cont, cont2, cont3, Nx, Ny)),:)=[cx cy cz];
           
        end
    end
end


for cont=1:2*(Nx+1)-1
    for cont2=1:2*(Ny+1)-1
        for cont3=1:Nz
            
            cx=0+sx*(ceil(cont/2)-1)+sx/2;
            cy=0+sy*(ceil(cont2/2)-1)+sy/2;
            cz=0+sz*(cont3-1)+sz/2;
            centri_vox.p1234((From_3D_to_1D(cont, cont2, cont3, 2*(Nx+1)-1, 2*(Ny+1)-1)),:)=[cx cy cz];
           
        end
    end
end


for cont=1:Nx
    for cont2=1:2*(Ny+1)-1
        for cont3=1:2*(Nz+1)-1
            
            cx=0+sx*(cont-1)+sx/2;
            cy=0+sy*(ceil(cont2/2)-1)+sy/2;
            cz=0+sz*(ceil(cont3/2)-1)+sz/2;
            centri_vox.p1256((From_3D_to_1D(cont, cont2, cont3, Nx, 2*(Ny+1)-1)),:)=[cx cy cz];
           
        end
    end
end

for cont=1:2*(Nx+1)-1
    for cont2=1:Ny
        for cont3=1:2*(Nz+1)-1
            
            cx=0+sx*(ceil(cont/2)-1)+sx/2;
            cy=0+sy*(cont2-1)+sy/2;
            cz=0+sz*(ceil(cont3/2)-1)+sz/2;
            centri_vox.p3456((From_3D_to_1D(cont, cont2, cont3, 2*(Nx+1)-1, Ny)),:)=[cx cy cz];
           
        end
    end
end

end