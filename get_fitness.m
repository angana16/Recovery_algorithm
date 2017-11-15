function [R,scores] = get_fitness(G,OD,indx,type)

reclist.edge_indx = indx;
reclist.nodes_indx = table2array(G.Edges(reclist.edge_indx,:));
reclist.nodes_indx(:,end) = [];

Gr = G;                             % original graph
Gr = rmedge(Gr,reclist.edge_indx);  % remove the edges in the rlist
scores = zeros(length(reclist.edge_indx),1); % initialize
sset = [];
for test = 1:length(reclist.edge_indx)
    scores(test) = inc_recover(Gr,sset,test,reclist.nodes_indx(test,1:2),OD,type);
    sset = [sset test];
    % update G
    Gr = addedge(Gr,reclist.nodes_indx(test,1),reclist.nodes_indx(test,2),1);
end

if length(scores) < length(reclist.edge_indx)
    temp1 = 1:length(reclist.edge_indx);
    temp2 = setdiff(temp1,sset);
    scores = [scores scores(end)*ones(1,length(temp2))];
    sset = [sset temp2];
end

nrlinks = length(scores);
R = nrlinks + 1 - ((scores(1) + scores(end))/2 + sum(scores(2:end-1)))/max(scores);
