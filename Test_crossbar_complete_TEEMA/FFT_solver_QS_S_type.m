function out=FFT_solver_QS_S_type(freq,escalings,incidence_selection,FFTCP,FFTCLp,...
    diagonals,ports,lumped_elements,expansions,GMRES_settings,Zs_info,QS_Rcc_FW)

warning off;

freq=freq*escalings.freq;

% GMRES settings ----------------------------
Inner_Iter=GMRES_settings.Inner_Iter;
Outer_Iter=GMRES_settings.Outer_Iter;
% -------------------------------------------

m=size(incidence_selection.A,1);
n=size(incidence_selection.A,2);
ns=size(incidence_selection.Gamma,2);

w=2*pi*freq;

nfreq=length(w);

is=zeros(n,1);

S=zeros(size(ports.port_nodes,1),size(ports.port_nodes,1),length(freq));

Vrest=zeros(m+n+ns,size(ports.port_nodes,1));

invP=sparse(1:ns,1:ns,1./diagonals.P,ns,ns);

R_chiusura=50;

for k=1:nfreq
    
    disp(['Freq n=' num2str(k) ' - Freq Tot=' num2str(nfreq)]);
    
    % Questa parte la vedremo in futuro (ignoratela per ora tanto QS_Rcc_FW è
    % impostato ad 1)-------------------------------------------------------
    if QS_Rcc_FW==2
        FFTCLp_rebuilted=compute_Circulant_Lp_Rcc(FFTCLp,escalings,freq(k)/escalings.freq);
        FFTCP_rebuilted=compute_Circulant_P_sup_Rcc(FFTCP,escalings,freq(k)/escalings.freq);
    elseif QS_Rcc_FW==3
        FFTCLp_rebuilted=compute_Circulant_Lp_FW(FFTCLp,escalings,freq(k)/escalings.freq);
        FFTCP_rebuilted=compute_Circulant_P_sup_FW(FFTCP,escalings,freq(k)/escalings.freq);
    end
    % ---------------------------------------------------------------------------
    
    Yle = build_Yle_S(lumped_elements,[],ports,escalings,n,w(k)/escalings.freq,R_chiusura);
    
    Z_self=compute_Z_self(diagonals.R,diagonals.Cd,w(k));
    
    
    Zs=escalings.R*(Zs_info.Zs*sqrt(w(k)/escalings.freq));
    
    ind_to_put_zero_Z_self=find((real(Zs(Zs_info.surface_edges))-real(Z_self(Zs_info.surface_edges)))>0);
    
    ind_to_put_zero_Zs=find((real(Zs(Zs_info.surface_edges))-real(Z_self(Zs_info.surface_edges)))<0);
    
    Z_self(Zs_info.surface_edges(ind_to_put_zero_Z_self))=0+1j*imag(Z_self(Zs_info.surface_edges(ind_to_put_zero_Z_self)));
    
    Zs(Zs_info.surface_edges(ind_to_put_zero_Zs))=0+1j*imag(Zs(Zs_info.surface_edges(ind_to_put_zero_Zs)));
    
    DZ=Z_self+real(Zs);
    
    DZ=DZ+1j*w(k)*diagonals.Lp;
    
    invZ=sparse(1:m,1:m,1./DZ,m,m);
    % --------------------- preconditioner ------------------------
    
    tic
    [L1,U1,P1,Q1]=lu((Yle+incidence_selection.A.'*invZ*incidence_selection.A+1j*w(k)*incidence_selection.Gamma*invP*incidence_selection.Gamma.'));
    
 
    time_Lu=toc;
    disp(['time LU ' num2str(time_Lu)]);
    
    % --------------------------------------------------------------
    
    
    for c1=1:size(ports.port_nodes,1)
        
        n1=ports.port_nodes(c1,1);
        n2=ports.port_nodes(c1,2);
        is(n1)=1*escalings.Is;
        is(n2)=-1*escalings.Is;
        
        tn=precond_3_3_Kt(L1,U1,P1,Q1,invZ,invP,incidence_selection,m,ns,is);
        
        
        if QS_Rcc_FW==1
            [V,flag,relres,iter,resvec]=...
                gmres(@ComputeMatrixVector,tn,Inner_Iter,GMRES_settings.tol(k),Outer_Iter,[],[],...
                Vrest(:,c1),w(k),incidence_selection,...
                FFTCP,FFTCLp,DZ,Yle,expansions,invZ,invP,L1,U1,P1,Q1);
            
        else
            [V,flag,relres,iter,resvec]=...
                gmres(@ComputeMatrixVector,tn,Inner_Iter,GMRES_settings.tol(k),Outer_Iter,[],[],...
                Vrest(:,c1),w(k),incidence_selection,...
                FFTCP_rebuilted,FFTCLp_rebuilted,DZ,Yle,expansions,invZ,invP,L1,U1,P1,Q1);
        end
        
        tot_iter_number=(iter(1)-1)*Inner_Iter+iter(2)+1;
        disp(['Flag ' num2str(flag) ' - Number of iterations = ' num2str(tot_iter_number) ]);
        
        Vrest(:,c1)=V;
        
        is(n1)=0;
        is(n2)=0;
        
        for c2=c1:size(ports.port_nodes,1)
            
            n3=ports.port_nodes(c2,1);
            n4=ports.port_nodes(c2,2);
            if c1==c2
                S(c1,c2,k)=(2*(V(m+ns+n3)-V(m+ns+n4))-R_chiusura)/R_chiusura;
            else
                S(c1,c2,k)=(2*(V(m+ns+n3)-V(m+ns+n4)))/R_chiusura;
            end
            S(c2,c1,k)=S(c1,c2,k);
            
        end
    end
    
end

out.S=S;
out.Z=s2z(S);
out.Y=s2y(S);
out.f=freq/escalings.freq;

end

function MatrixVector=ComputeMatrixVector(x,w,incidence_selection,...
    FFTCP,FFTCLp,DZ,Yle,expansions,invZ,invP,L1,U1,P1,Q1)

m=size(incidence_selection.A,1);
ns=size(incidence_selection.Gamma,2);

I=x(1:m);
Q=x(m+1:m+ns);
Phi=x(m+ns+1:end);

% Lp * I ---------------------------------------------------------------
mx=incidence_selection.mx;
my=incidence_selection.my;
mz=incidence_selection.mz;

Y1=zeros(m,1);

ind_aux_Lp{1}=1:mx;
ind_aux_Lp{2}=mx+1:mx+my;
ind_aux_Lp{3}=mx+my+1:mx+my+mz;

for cont=1:3
    
    Nx=size(FFTCLp{cont,1},1)/2;
    Ny=size(FFTCLp{cont,1},2)/2;
    Nz=size(FFTCLp{cont,1},3)/2;
    Ired=I(ind_aux_Lp{cont});
    I_exp=expansions.mat_map_Lp{cont,1}*Ired;
    CircKT=reshape(I_exp,Nx,Ny,Nz);
    Chi=ifftn(FFTCLp{cont,1}.*fftn(CircKT,[2*Nx,2*Ny,2*Nz]));
    Y1(ind_aux_Lp{cont})=Y1(ind_aux_Lp{cont})+expansions.mat_map_Lp{cont,1}.'*reshape(Chi(1:Nx,1:Ny,1:Nz),Nx*Ny*Nz,1);
    
end

Y1=1j*w*Y1+DZ.*I+incidence_selection.A*Phi;


% ---------------- P * Q ---------------------------------------------

Y2=zeros(ns,1);

for cont1=1:3
    for cont2=cont1:3
        Nx=size(FFTCP{cont1,cont2},1)/2;
        Ny=size(FFTCP{cont1,cont2},2)/2;
        Nz=size(FFTCP{cont1,cont2},3)/2;
        
        Q_exp=expansions.exp_P{cont1,cont2}*Q;
        CircKT=reshape(Q_exp,Nx,Ny,Nz);
        Chi=ifftn(FFTCP{cont1,cont2}.*fftn(CircKT,[2*Nx,2*Ny,2*Nz]));
        Y2=Y2+expansions.exp_P{cont2,cont1}.'*reshape(Chi(1:Nx,1:Ny,1:Nz),Nx*Ny*Nz,1);
        
        if cont1~=cont2
            Q_exp=expansions.exp_P{cont2,cont1}*Q;
            CircKT=reshape(Q_exp,Nx,Ny,Nz);
            Chi=ifftn(FFTCP{cont1,cont2}.*fftn(CircKT,[2*Nx,2*Ny,2*Nz]));
            Y2=Y2+expansions.exp_P{cont1,cont2}.'*reshape(Chi(1:Nx,1:Ny,1:Nz),Nx*Ny*Nz,1);
        end
    end
end

Y2=Y2-incidence_selection.Gamma.'*Phi;

Y3=-incidence_selection.A.'*I+Yle*Phi+1j*w*(incidence_selection.Gamma*Q);

MatrixVector=precond_3_3_vector(L1,U1,P1,Q1,invZ,invP,incidence_selection,w,Y1,Y2,Y3);

end

function Y=precond_3_3_Kt(L1,U1,P1,Q1,invZ,invP,incidence_selection,n1,n2,X3)

n3=length(X3);

i1=1:n1;
i2=n1+1:n1+n2;
i3=n1+n2+1:n1+n2+n3;

Y=zeros(n1+n2+n3,1);

M5 = Q1*(U1\(L1\(P1*X3)));

Y(i1)=Y(i1)-invZ*(incidence_selection.A*M5);

Y(i2)=Y(i2)+invP*(incidence_selection.Gamma.'*M5);

Y(i3)=Y(i3)+M5;

end

function Y=precond_3_3_vector(L1,U1,P1,Q1,invZ,invP,incidence_selection,w,X1,X2,X3)

n1=length(X1);
n2=length(X2);
n3=length(X3);

i1=1:n1;
i2=n1+1:n1+n2;
i3=n1+n2+1:n1+n2+n3;

Y=zeros(n1+n2+n3,1);

M1 = invZ*X1;
M2 = (Q1*(U1\(L1\(P1*(incidence_selection.A.'*M1)))));
M3 = invP*X2;
M4 = (Q1*(U1\(L1\(P1*(incidence_selection.Gamma*M3)))));
M5 = Q1*(U1\(L1\(P1*X3)));

Y(i1)=Y(i1)+invZ*X1-invZ*(incidence_selection.A*M2);
Y(i1)=Y(i1)+1j*w*(invZ*(incidence_selection.A*M4));
Y(i1)=Y(i1)-invZ*(incidence_selection.A*M5);

Y(i2)=Y(i2)+invP*(incidence_selection.Gamma.'*M2);
Y(i2)=Y(i2)+invP*X2-1j*w*invP*(incidence_selection.Gamma.'*M4);
Y(i2)=Y(i2)+invP*(incidence_selection.Gamma.'*M5);

Y(i3)=Y(i3)+M2;
Y(i3)=Y(i3)-1j*w*M4;
Y(i3)=Y(i3)+M5;

end


