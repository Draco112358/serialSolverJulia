using Distances

function find_voxels_port_pp(nodi_centri, ports, nodes, nodes_red)
    N = size(ports.port_start, 1)
    ports.port_nodes = zeros(N, 2)
    ports.port_voxels = zeros(N, 2)
    for cont = 1:N
        ports.port_voxels[cont, 1] = nodes_find_rev(ports.port_start[cont, :], nodi_centri, -1)
        ports.port_voxels[cont, 2] = nodes_find_rev(ports.port_end[cont, :], nodi_centri, ports.port_voxels[cont, 1])
        ports.port_nodes[cont, 1] = searchsorted(nodes[ports.port_voxels[cont, 1]], nodes_red)
        ports.port_nodes[cont, 2] = searchsorted(nodes[ports.port_voxels[cont, 2]], nodes_red)
    end
end

function nodes_find_rev(Nodes_inp_coord, nodi_centri, node_to_skip)
    indici = sort(collect(distances(Nodes_inp_coord', nodi_centri')))[2:end]
    if indici[1] != node_to_skip
        nodes = indici[1]
    else
        nodes = indici[2]
    end
end
