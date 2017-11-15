function p = plotting(adj,nodes)
G = graph(adj);
G.Nodes.id = nodes.ID;
nnodes = size(G.Nodes,1);
nlinks = size(G.Edges,1);

% nodes_snow = nodes_snow2015;
% plot network with the damaged nodes
%------------------------------------------------
f = figure;
p = plot(G,'layout','force');
p.Marker = 'o';
p.NodeColor = [0 0 0];
% p.NodeFaceColor = [1 1 1];
p.EdgeColor = [0 0 0];
%temp = 5*ones(nnodes,1);
%temp(nodes_snow.index) = 5;
%G.Nodes.Snow = temp;
%p.MarkerSize = temp;
p.XData = nodes.LONGITUDE;
p.YData = nodes.LATITUDE;
end
