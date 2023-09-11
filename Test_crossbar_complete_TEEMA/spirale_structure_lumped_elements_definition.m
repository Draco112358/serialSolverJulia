function lumped_elements=spirale_structure_lumped_elements_definition()

% define the lumped elements position in meters
num_le=1;
le_end=zeros(num_le,3);
le_start=zeros(num_le,3);
lumped_elements.value=[1e12];
lumped_elements.type=ones(num_le,1);

le_start(1,:)=[-1.9 -0.344 0.025]*1e-3;
le_end(1,:)=[-1.9 -0.5965 0.025]*1e-3;

lumped_elements.le_start=le_start;
lumped_elements.le_end=le_end;

lumped_elements.s_le_start=[];
lumped_elements.s_le_end=[];

end