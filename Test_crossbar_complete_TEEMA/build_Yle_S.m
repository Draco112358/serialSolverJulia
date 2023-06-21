function Yle = build_Yle_S(lumped_elements,grounding_nodes,ports,escalings,n,w,val_chiusura)

% da fare per i vari tipi di circuito

le_nodes=[lumped_elements.le_nodes(:,1);lumped_elements.le_nodes(:,2)];
port_nodes=[ports.port_nodes(:,1);ports.port_nodes(:,2)];
N_ele=length(unique([le_nodes;port_nodes;grounding_nodes.']));

NNz_max=N_ele^2;

ind_r=zeros(NNz_max,1);
ind_c=zeros(NNz_max,1);
vals=zeros(NNz_max,1);

nlum=size(lumped_elements.le_nodes);

cont=0;
for c1=1:nlum
    
    n1=lumped_elements.le_nodes(c1,1);
    n2=lumped_elements.le_nodes(c1,2);
    
    ind1=find(ind_r==n1);
    ind2=find(ind_c==n1);
    
    ind=intersect(ind1,ind2);
    
    if isempty(ind)
        cont=cont+1;
        ind_r(cont)=n1;
        ind_c(cont)=n1;
        if lumped_elements.type(c1)==3
            vals(cont)=1/(1j*w*lumped_elements.value(c1));
        elseif lumped_elements.type(c1)==2
            vals(cont)=1j*w*lumped_elements.value(c1);
        else
            vals(cont)=1/lumped_elements.value(c1);
        end
    else
        if lumped_elements.type(c1)==3
            vals(ind(1))=vals(ind(1))+1/(1j*w*lumped_elements.value(c1));
        elseif lumped_elements.type(c1)==2
            vals(ind(1))=vals(ind(1))+1j*w*lumped_elements.value(c1);
        else
            vals(ind(1))=vals(ind(1))+1/lumped_elements.value(c1);
        end
    end
    
    ind1=find(ind_r==n2);
    ind2=find(ind_c==n2);
    
    ind=intersect(ind1,ind2);
    
    if isempty(ind)
        cont=cont+1;
        ind_r(cont)=n2;
        ind_c(cont)=n2;
        if lumped_elements.type(c1)==3
            vals(cont)=1/(1j*w*lumped_elements.value(c1));
        elseif lumped_elements.type(c1)==2
            vals(cont)=1j*w*lumped_elements.value(c1);
        else
            vals(cont)=1/lumped_elements.value(c1);
        end
    else
        if lumped_elements.type(c1)==3
            vals(ind(1))=vals(ind(1))+1/(1j*w*lumped_elements.value(c1));
        elseif lumped_elements.type(c1)==2
            vals(ind(1))=vals(ind(1))+1j*w*lumped_elements.value(c1);
        else
            vals(ind(1))=vals(ind(1))+1/lumped_elements.value(c1);
        end
    end
    
    ind1=find(ind_r==n1);
    ind2=find(ind_c==n2);
    
    ind=intersect(ind1,ind2);
    
    if isempty(ind)
        cont=cont+1;
        ind_r(cont)=n1;
        ind_c(cont)=n2;
        if lumped_elements.type(c1)==3
            vals(cont)=-1/(1j*w*lumped_elements.value(c1));
        elseif lumped_elements.type(c1)==2
            vals(cont)=-1j*w*lumped_elements.value(c1);
        else
            vals(cont)=-1/lumped_elements.value(c1);
        end
    else
        if lumped_elements.type(c1)==3
            vals(ind(1))=vals(ind(1))-1/(1j*w*lumped_elements.value(c1));
        elseif lumped_elements.type(c1)==2
            vals(ind(1))=vals(ind(1))-1j*w*lumped_elements.value(c1);
        else
            vals(ind(1))=vals(ind(1))-1/lumped_elements.value(c1);
        end
    end
    
    ind1=find(ind_r==n2);
    ind2=find(ind_c==n1);
    
    ind=intersect(ind1,ind2);
    
    if isempty(ind)
        cont=cont+1;
        ind_r(cont)=n2;
        ind_c(cont)=n1;
        if lumped_elements.type(c1)==3
            vals(cont)=-1/(1j*w*lumped_elements.value(c1));
        elseif lumped_elements.type(c1)==2
            vals(cont)=-1j*w*lumped_elements.value(c1);
        else
            vals(cont)=-1/lumped_elements.value(c1);
        end
    else
        if lumped_elements.type(c1)==3
            vals(ind(1))=vals(ind(1))-1/(1j*w*lumped_elements.value(c1));
        elseif lumped_elements.type(c1)==2
            vals(ind(1))=vals(ind(1))-1j*w*lumped_elements.value(c1);
        else
            vals(ind(1))=vals(ind(1))-1/lumped_elements.value(c1);
        end
    end
    
end

nlum=size(ports.port_nodes);

for c1=1:nlum
    
    n1=ports.port_nodes(c1,1);
    n2=ports.port_nodes(c1,2);
    
    ind1=find(ind_r==n1);
    ind2=find(ind_c==n1);
    
    ind=intersect(ind1,ind2);
    
    if isempty(ind)
        cont=cont+1;
        ind_r(cont)=n1;
        ind_c(cont)=n1;
        vals(cont)=1/val_chiusura;
    else
        vals(ind(1))=vals(ind(1))+1/val_chiusura;
    end
    
    ind1=find(ind_r==n2);
    ind2=find(ind_c==n2);
    
    ind=intersect(ind1,ind2);
    
    if isempty(ind)
        cont=cont+1;
        ind_r(cont)=n2;
        ind_c(cont)=n2;
         vals(cont)=1/val_chiusura;
    else
        vals(ind(1))=vals(ind(1))+1/val_chiusura;
    end
    
    ind1=find(ind_r==n1);
    ind2=find(ind_c==n2);
    
    ind=intersect(ind1,ind2);
    
    if isempty(ind)
        cont=cont+1;
        ind_r(cont)=n1;
        ind_c(cont)=n2;
         vals(cont)=-1/val_chiusura;
    else
        vals(ind(1))=vals(ind(1))-1/val_chiusura;
    end
    
    ind1=find(ind_r==n2);
    ind2=find(ind_c==n1);
    
    ind=intersect(ind1,ind2);
    
    if isempty(ind)
        cont=cont+1;
        ind_r(cont)=n2;
        ind_c(cont)=n1;
         vals(cont)=-1/val_chiusura;
    else
        vals(ind(1))=vals(ind(1))-1/val_chiusura;
    end
    
end

for k=1:length(grounding_nodes)
    
    n1=grounding_nodes(k);
    
    ind1=find(ind_r==n1);
    ind2=find(ind_c==n1);
    
    ind=intersect(ind1,ind2);
    
    if isempty(ind)
        cont=cont+1;
        vals(cont)=1/1e12;
        ind_r(cont)=grounding_nodes(k);
        ind_c(cont)=grounding_nodes(k);
    else
        vals(ind(1))=vals(ind(1))+1/1e12;
    end
    
end

Yle = sparse(ind_r(1:cont),ind_c(1:cont),vals(1:cont)*escalings.Yle,n,n);
end