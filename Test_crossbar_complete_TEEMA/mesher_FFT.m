function [escalings,incidence_selection,circulant_centers,diagonals,expansions,ports,...
    lumped_elements,li_mats,Zs_info]=mesher_FFT(use_escalings,materials,sx,sy,sz,grids,centri_vox,externals_grids,mapping_vols,ports,lumped_elements)

global min_v dominant_list

if isempty(dominant_list)
    dominant_list=1;
end

escalings.Lp=1;
escalings.P=1;
escalings.R=1;
escalings.Cd=1;
escalings.Is=1;
escalings.Yle=1;
escalings.freq=1;
escalings.time=1;

if use_escalings==1
    
    escalings.Lp=1e6;
    escalings.P=1e-12;
    escalings.R=1e-3;
    escalings.Cd=1e12;
    escalings.Is=1e3;
    escalings.Yle=1e3;
    escalings.freq=1e-9;
    escalings.time=1e9;
    
end

Nx=size(grids{1,1},1);
Ny=size(grids{1,1},2);
Nz=size(grids{1,1},3);

num_input_files=size(materials,1);

for cont=1:num_input_files
    materials{cont}.epsr=materials{cont}.eps_re-1j*materials{cont}.eps_re*materials{cont}.tan_D;
end

num_full_vox=0;
for cont=1:num_input_files
    num_full_vox=num_full_vox+nnz(grids{cont});
end

[nodes,nodes_red,nodes_reused_clean]=create_nodes_ref(grids,num_full_vox,externals_grids,mapping_vols,dominant_list);

[expansions,incidence_selection.Gamma]=create_mapping_Gamma_no_rep(grids,mapping_vols,nodes,nodes_red,externals_grids);

disp(['Total Full surfaces          ' num2str(size(incidence_selection.Gamma,2))]);

[mappings.Ax,mappings.NAx]=create_mapping_Ax_v2(grids,mapping_vols,nodes,nodes_red);
[mappings.Ay,mappings.NAy]=create_mapping_Ay_v2(grids,mapping_vols,nodes,nodes_red);
[mappings.Az,mappings.NAz]=create_mapping_Az_v2(grids,mapping_vols,nodes,nodes_red);

[incidence_selection.A,li_mats.lix_mat,li_mats.liy_mat,li_mats.liz_mat,li_mats.lix_border,li_mats.liy_border,li_mats.liz_border,maps_Zs,...
    ]=create_A_mats_and_find_borders_with_map_Zs(grids,mapping_vols,...
    mappings.Ax, mappings.NAx,mappings.Ay, mappings.NAy,mappings.Az, mappings.NAz, materials,nodes,nodes_red,nodes_reused_clean);

[Zs_info,ind_cond_rug_x,ind_cond_rug_y,ind_cond_rug_z]=compute_Zs_with_indices(materials,li_mats,maps_Zs,sx,sy,sz);
Zs_info.ind_cond_rug_x=ind_cond_rug_x;
Zs_info.ind_cond_rug_y=ind_cond_rug_y;
Zs_info.ind_cond_rug_z=ind_cond_rug_z;
clear ind_cond_rug_x ind_cond_rug_y ind_cond_rug_z

disp(['Total inductive edges        ' num2str(size(incidence_selection.A,1))]);
disp(['Total  nodes                 ' num2str(size(incidence_selection.A,2))]);

incidence_selection.mx=size(li_mats.lix_mat,1);
incidence_selection.my=size(li_mats.liy_mat,1);
incidence_selection.mz=size(li_mats.liz_mat,1);

ind_Lp=cell(3,1);
expansions.N_ind=cell(3,1);
expansions.mat_map_Lp=cell(3,3);

ind_Lp{1}=create_expansion_ind_Lp_x_grids_v2(grids,mappings.Ax,mappings.NAx,mapping_vols,nodes,nodes_red);
expansions.N_ind{1}=(Nx-1)*Ny*Nz;

expansions.mat_map_Lp{1,1}=sparse(ind_Lp{1,1}(:,1),...
    ind_Lp{1,1}(:,2),...
    ones(size(ind_Lp{1,1},1),1),...
    expansions.N_ind{1,1},size(ind_Lp{1,1},1));


lati_tenere=li_mats.lix_border(:,1)~=0;
ind_aux=ind_Lp{1,1}(lati_tenere,:);

expansions.mat_map_Lp{1,2}=sparse(ind_aux(:,1),...
    ind_aux(:,2),...
    ones(size(ind_aux,1),1),...
    expansions.N_ind{1,1},size(ind_Lp{1},1));

lati_tenere=li_mats.lix_border(:,2)~=0;
ind_aux=ind_Lp{1,1}(lati_tenere,:);

expansions.mat_map_Lp{1,3}=sparse(ind_aux(:,1),...
    ind_aux(:,2),...
    ones(size(ind_aux,1),1),...
    expansions.N_ind{1,1},size(ind_Lp{1},1));

ind_Lp{2}=create_expansion_ind_Lp_y_grids_v2(grids,mappings.Ay,mappings.NAy,mapping_vols,nodes,nodes_red);
expansions.N_ind{2}=(Ny-1)*Nx*Nz;

expansions.mat_map_Lp{2,1}=sparse(ind_Lp{2}(:,1),...
    ind_Lp{2}(:,2),...
    ones(size(ind_Lp{2},1),1),...
    expansions.N_ind{2,1},size(ind_Lp{2},1));

lati_tenere=li_mats.liy_border(:,1)~=0;
ind_aux=ind_Lp{2}(lati_tenere,:);

expansions.mat_map_Lp{2,2}=sparse(ind_aux(:,1),...
    ind_aux(:,2),...
    ones(size(ind_aux,1),1),...
    expansions.N_ind{2,1},size(ind_Lp{2},1));

lati_tenere=li_mats.liy_border(:,2)~=0;
ind_aux=ind_Lp{2}(lati_tenere,:);

expansions.mat_map_Lp{2,3}=sparse(ind_aux(:,1),...
    ind_aux(:,2),...
    ones(size(ind_aux,1),1),...
    expansions.N_ind{2,1},size(ind_Lp{2},1));

ind_Lp{3}=create_expansion_ind_Lp_z_grids_v2(grids,mappings.Az,mappings.NAz,mapping_vols,nodes,nodes_red);
expansions.N_ind{3}=(Nz-1)*Nx*Ny;

expansions.mat_map_Lp{3,1}=sparse(ind_Lp{3}(:,1),...
    ind_Lp{3}(:,2),...
    ones(size(ind_Lp{3},1),1),...
    expansions.N_ind{3,1},size(ind_Lp{3},1));

lati_tenere=li_mats.liz_border(:,1)~=0;
ind_aux=ind_Lp{3}(lati_tenere,:);

expansions.mat_map_Lp{3,2}=sparse(ind_aux(:,1),...
    ind_aux(:,2),...
    ones(size(ind_aux,1),1),...
    expansions.N_ind{3,1},size(ind_Lp{3},1));

lati_tenere=li_mats.liz_border(:,2)~=0;
ind_aux=ind_Lp{3}(lati_tenere,:);

expansions.mat_map_Lp{3,3}=sparse(ind_aux(:,1),...
    ind_aux(:,2),...
    ones(size(ind_aux,1),1),...
    expansions.N_ind{3,1},size(ind_Lp{3},1));

diagonals=compute_diagonals(escalings,materials,sx,sy,sz,li_mats.lix_mat,li_mats.liy_mat,li_mats.liz_mat,li_mats.lix_border,li_mats.liy_border,li_mats.liz_border);

li_mats.sx=sx;
li_mats.sy=sy;
li_mats.sz=sz;
diagonals.P=compute_diagonal_P_v2(expansions.N1,expansions.N2,expansions.N3,escalings,sx,sy,sz);


% ---only point-to-point ports or le ------------------
ports=find_voxels_port_pp(centri_vox,ports,nodes,nodes_red);
ports.surf_s_port_nodes=cell(0,1);
ports.surf_e_port_nodes=cell(0,1);
lumped_elements=find_voxels_le_pp(centri_vox,lumped_elements,nodes,nodes_red);
lumped_elements.surf_s_le_nodes=zeros(0,3);
lumped_elements.surf_e_le_nodes=zeros(0,3);
% ---------------------

circulant_centers=build_center_P_Voxels(min_v,Nx,Ny,Nz,sx,sy,sz);

[circulant_centers.Lpx,circulant_centers.Lpy,circulant_centers.Lpz]=build_centers_Lp_with_air(grids,build_center_air_Voxels(grids,sx,sy,sz,min_v));

circulant_centers.sx=sx;
circulant_centers.sy=sy;
circulant_centers.sz=sz;

circulant_centers.Nx=Nx;
circulant_centers.Ny=Ny;
circulant_centers.Nz=Nz;

end