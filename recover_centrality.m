function [scores,reclist] = recover_centrality(G,OD,ucc,type,rnodes_indx,list)

temp1 = ucc(rnodes_indx);
[~,indx] = sort(temp1,'descend');
sset_nodes = rnodes_indx(indx);        % order of recovery nodes

indx = [];
for i = 1:length(sset_nodes)
    var2 = ismember(list.nodes_indx,sset_nodes(i));
    var3 = +(sum(var2,2)>0);
    indx = [indx; find(var3)];
end
[indx] = unique(indx,'rows','stable');  % remove duplicates
reclist.edge_indx = indx;
reclist.nodes_indx = table2array(G.Edges(reclist.edge_indx,:));
reclist.nodes_indx(:,end) = [];

Gr = G;                             % original graph
Gr = rmedge(Gr,reclist.edge_indx);  % remove the edges in the rlist
scores = zeros(length(reclist.edge_indx),1); % initialize
sset = [];
reclist
for test = 1:length(reclist.edge_indx)
    scores(test) = inc_recover(Gr,sset,test,reclist.nodes_indx(test,1:2),OD,type);
    sset = [sset test];
    % update G
    Gr = addedge(Gr,reclist.nodes_indx(test,1),reclist.nodes_indx(test,2),1);
end

