function nodes_ports=crossbar_250_structure_nodes_ports_definition()

% define the port position in meters

num_ports=2;
port_start=zeros(num_ports,3);
port_end=zeros(num_ports,3);

k=1;
port_start(k,:)=[0 0.375 0.0312]*1e-3;
port_end(k,:)=  [0 0.475 0.0312]*1e-3;k=k+1;

port_start(k,:)=[0.375 0 0.0063]*1e-3;
port_end(k,:)=  [0.475 0 0.0063]*1e-3;k=k+1;

nodes_ports.port_start=port_start;
nodes_ports.port_end=port_end;

nodes_ports.s_port_start=zeros(0,6);
nodes_ports.s_port_end=zeros(0,6);

end