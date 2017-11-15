function [adjacency,weighted_adjacency,node_struc] = adjacency_create(nodes, edges)
% run as [adjacency,weighted_adjacency,node_struc] = adjacency_create(US_nodes, US_edges)
node_struc.index = nodes.ID;
node_struc.id = nodes.Label;
node_struc.lat = nodes.LATITUDE;
node_struc.lon = nodes.LONGITUDE;
edge_struc.source = edges.Source;
edge_struc.target = edges.Target;
edge_struc.weight = edges.Weight;

adjacency =  zeros(length(node_struc.id),length(node_struc.id));
weighted_adjacency = zeros(length(node_struc.id),length(node_struc.id));

for i =1:length(edge_struc.source)
adjacency(edge_struc.source(i),edge_struc.target(i)) = 1;
adjacency(edge_struc.target(i),edge_struc.source(i)) = 1;

weighted_adjacency(edge_struc.source(i),edge_struc.target(i)) = edge_struc.weight(i);
weighted_adjacency(edge_struc.target(i),edge_struc.source(i)) = edge_struc.weight(i);


end 
end 
