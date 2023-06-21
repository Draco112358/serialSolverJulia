function Z_self=compute_Z_self(R,Cd,w)

Z_self=zeros(length(R),1)+1j*zeros(length(R),1);
for cont=1:length(R)
    
    for aux=1:4
        
        if R(cont,aux)~=0 && Cd(cont,aux)~=0
            Z_self(cont)=Z_self(cont)+1/(1/R(cont,aux)+1j*w*Cd(cont,aux));
        elseif R(cont,aux)~=0
            Z_self(cont)=Z_self(cont)+R(cont,aux);
        elseif Cd(cont,aux)~=0
            Z_self(cont)=Z_self(cont)+1/(1j*w*Cd(cont,aux));
        end
    end
end

end