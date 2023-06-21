function [diag_P]=compute_diagonal_P_v2(N1,N2,N3,escalings,sx,sy,sz)
 
self_P=zeros(3,1);
centro_vox=[0 0 0];
row_P=escalings.P*compute_row_P_sup(centro_vox,centro_vox,sx,sy,sz,1,1);
self_P(1)=row_P(1);
row_P=escalings.P*compute_row_P_sup(centro_vox,centro_vox,sx,sy,sz,3,3);
self_P(2)=row_P(1);
row_P=escalings.P*compute_row_P_sup(centro_vox,centro_vox,sx,sy,sz,5,5);
self_P(3)=row_P(1);

diag_P=[...
    self_P(1)*ones(N1,1);...
    self_P(2)*ones(N2,1);...
    self_P(3)*ones(N3,1);...
];

end