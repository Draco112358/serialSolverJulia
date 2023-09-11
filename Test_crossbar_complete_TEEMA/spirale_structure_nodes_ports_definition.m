function nodes_ports=spirale_structure_nodes_ports_definition()

% define the port position in meters

num_ports=1;
port_start=zeros(num_ports,3);
port_end=zeros(num_ports,3);

k=1;
port_start(k,:)=[-1.9 -0.344 0.025]*1e-3;
port_end(k,:)=  [-1.9 -0.5965 0.025]*1e-3;

nodes_ports.port_start=port_start;
nodes_ports.port_end=port_end;

nodes_ports.s_port_start=zeros(0,6);
nodes_ports.s_port_end=zeros(0,6);

end