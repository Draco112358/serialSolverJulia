clc
close all
clearvars
clear global

% input setup -------------------------------------------------------------
% k=1;
% materials{k}.sigmar=9.4e6;   materials{k}.eps_re=1;      materials{k}.tan_D=0; materials{k}.mur=1;k=k+1;
% 
% ports=crossbar_structure_nodes_ports_definition();
% lumped_elements=crossbar_structure_lumped_elements_definition();

%ports=spirale_structure_nodes_ports_definition();
%lumped_elements=spirale_structure_lumped_elements_definition();
global minmax

%load('mesher_output.mat');
%load('/Users/edgardovittoria/Downloads/mesher_output_spirale.mat')
load('/Users/edgardovittoria/Downloads/matfile.mat')
% jsonstring = jsonencode(grids);
% fid = fopen('/tmp/mehseroutputmatlab.json','w');
% fprintf(fid,'%s',jsonstring);
% fclose(fid);

% x = numel(grids{1})
% y = numel(grids{1}{1})
% z = numel(grids{1}{1}{1})
% 
% multiDimCell = true(x, y, z);
% for i=1:x
%     for j=1:y
%         for k=1:z
%             disp(k)
%             multiDimCell(i,j,k) = grids{1,1}{i}{j}(k)
%         end
%     end
% end
% 
% grids = cell(1,1)
% grids{1} = multiDimCell
% 
% save('/Users/edgardovittoria/Downloads/mesher_output_spirale.mat', "grids", "num_full_vox", "sx", "sy", "sz")

% Nx=size(grids{1},1);
% Ny=size(grids{1},2);
% Nz=size(grids{1},3);
% 
% % -------------------------------------------------------------------------
% 
% % -------Solver setup ---------------------------------------------------
% n_freq=10;
% freq=logspace(log10(1e6),log10(1e9),n_freq);
% 
% GMRES_settings.Inner_Iter=150;
% GMRES_settings.Outer_Iter=1;
% 
% GMRES_settings.tol=1e-4*ones(n_freq,1);
% ind_low_freq=find(freq<1e5);
% GMRES_settings.tol(ind_low_freq)=1e-4;
% 
% QS_Rcc_FW=1; % 1 QS, 2 Rcc, 3 Taylor
% use_escalings=1;
% % -------------------------------------------------------------------------
% [mapping_vols,num_centri]=create_volumes_mapping_v2(grids);
% [centri_vox,id_mat]=create_volume_centers(grids,mapping_vols,num_full_vox,sx,sy,sz);
% 
% externals_grids=create_Grids_externals(grids);
% 
% [escalings,incidence_selection,circulant_centers,diagonals,expansions,ports,...
%     lumped_elements,li_mats,Zs_info]=mesher_FFT(use_escalings,...
%     materials,sx,sy,sz,grids,centri_vox,externals_grids,mapping_vols,ports,lumped_elements);
% 
% [FFTCP,FFTCLp]=compute_FFT_mutual_coupling_mats(circulant_centers,escalings,Nx,Ny,Nz,QS_Rcc_FW);

 out=FFT_solver_QS_S_type(freq,escalings,incidence_selection,FFTCP,FFTCLp,...
    diagonals,ports,lumped_elements,expansions,GMRES_settings,Zs_info,QS_Rcc_FW);
 
figure
semilogx(freq,real(squeeze(out.Z(1,1,:))),'b*',...
    'linewidth',2)
xlim([freq(1) 1e9])
xlabel('Frequency [Hz]','fontsize',14)
ylabel(['R [\Omega]'],'fontsize',14)
legend('MATLAB')
h_legend=legend;
set(gca,'XTick',[10 10^2 10^3 10^4 10^5 10^6 10^7 10^8 10^9])
xticklabels({'10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}','10^{7}','10^{8}','10^{9}'})
set(h_legend,'Location','northwest','FontSize',14);
set(gca,'fontsize',14);

figure
semilogx(freq,imag(squeeze(out.Z(1,1,:)))./(2*pi*freq.')*1e9,'b*',...
    'linewidth',2)
xlim([freq(1) 1e9])
xlabel('Frequency [Hz]','fontsize',14)
ylabel(['L [nH]'],'fontsize',14)
legend('MATLAB')
h_legend=legend;
set(gca,'XTick',[10 10^2 10^3 10^4 10^5 10^6 10^7 10^8 10^9])
xticklabels({'10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}','10^{7}','10^{8}','10^{9}'})
set(h_legend,'Location','northwest','FontSize',14);
set(gca,'fontsize',14);


figure
semilogx(freq,real(squeeze(out.Z(1,2,:))),'b*',...
    'linewidth',2)
xlim([freq(1) 1e9])
xlabel('Frequency [Hz]','fontsize',14)
ylabel(['R [\Omega]'],'fontsize',14)
legend('MATLAB')
h_legend=legend;
set(gca,'XTick',[10 10^2 10^3 10^4 10^5 10^6 10^7 10^8 10^9])
xticklabels({'10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}','10^{7}','10^{8}','10^{9}'})
set(h_legend,'Location','northwest','FontSize',14);
set(gca,'fontsize',14);

figure
semilogx(freq,imag(squeeze(out.Z(1,2,:)))./(2*pi*freq.')*1e9,'b*',...
    'linewidth',2)
xlim([freq(1) 1e9])
xlabel('Frequency [Hz]','fontsize',14)
ylabel(['L [nH]'],'fontsize',14)
legend('MATLAB')
h_legend=legend;
set(gca,'XTick',[10 10^2 10^3 10^4 10^5 10^6 10^7 10^8 10^9])
xticklabels({'10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}','10^{7}','10^{8}','10^{9}'})
set(h_legend,'Location','northwest','FontSize',14);
set(gca,'fontsize',14);


