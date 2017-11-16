% This is master code to generate recovery profiles. No need to modify any
% other function. Make sure all other functions are added in your working
% directory. 
function [c,scores,T] = resilience_CIR(Adj,adj_weighted, nodes,rm_nodes,type )
% Written by: Udit Bhatia and Lina Sela
% Modified by: Udit Bhatia
% Date last modified: November 13, 2017
% Run as following for indian railways in command window (load Indian_datafiles.mat):
%[c,scores,T ] = resilience_CIR(adj_IRN,adj_IRN_weighted, nodes_irn,tsunami,"LargeC")
%Run as following for MBTA in command window (load mbta_data.mat):
%[c,scores,T ] = resilience_CIR(Adj,Adj_weighted, nodes,nodes_snow,"LargeC")
%Run as following for USNAS in command window (load US_datafiles.mat):
%[c,scores,T ] = resilience_CIR(adjacency,weighted_adjacency, node_struc,nodes_rem,"LargeC")

G = graph(Adj);
G.Nodes.id = nodes.id;
nnodes = size(G.Nodes,1);
nnlinks = size (G.Nodes,1);
var1 = rm_nodes.index;
list.nodes_indx = table2array(G.Edges);     
list.nodes_indx(:,end) = [];               
var2 = ismember(list.nodes_indx,var1);

var3 = +(sum(var2,2)>0);
rlist.edge_indx = find(var3);
rlist.nodes_indx = table2array(G.Edges(rlist.edge_indx,:));
rlist.nodes_indx(:,end) = [];
budget = length(rlist.edge_indx);  % Number of edges removed

OD = adj_weighted;
[greedy.sset,greedy.scores,greedy.evalNum] = greedy_lazy(G, OD, rlist, budget,type);
Gr = graph(Adj);
scores = zeros(length(rlist.edge_indx),1);
for i = 1:length(rlist.edge_indx)
    Gr = rmedge(Gr,rlist.nodes_indx(i,1),rlist.nodes_indx(i,2));
    scores(i) = ODScore(Gr,OD,type);
end
if length(greedy.scores) < length(rlist.edge_indx)
    temp1 = 1:length(rlist.edge_indx);
    temp2 = setdiff(temp1,greedy.sset);
    greedy.scores = [greedy.scores greedy.scores(end)*ones(1,length(temp2))];
    greedy.sset = [greedy.sset temp2];
end
greedy.sset = rlist.edge_indx(greedy.sset);

Scores = [scores; greedy.scores'];
greedy_scores = Scores;

c.ucc = centrality(G,'closeness');
c.ud = centrality(G,'degree');
c.ubc = centrality(G,'betweenness');
c.ubc = 2*c.ubc/((nnodes-2)*(nnodes-1));
c.ueig = centrality(G,'eigenvector');

rnodes_indx = rm_nodes.index;        % removed nodes index
% closenesspip in
[c.scores.ucc,c.reclist.ucc] = recover_centrality(G,OD,c.ucc,type,rnodes_indx,list);
[c.scores.ud,c.reclist.ud] = recover_centrality(G,OD,c.ud,type,rnodes_indx,list);
[c.scores.ubc,c.reclist.ubc] = recover_centrality(G,OD,c.ubc,type,rnodes_indx,list);
[c.scores.ueig,c.reclist.ueig] = recover_centrality(G,OD,c.ueig,type,rnodes_indx,list);

f = figure();

Scores2 = [scores; c.scores.ucc]; % combine failure and recovery scores
hold on
plot(Scores2/max(Scores2),'m','linewidth',3)
Scores2 = [scores; c.scores.ud]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'k','linewidth',3)
Scores2 = [scores; c.scores.ubc]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'g','linewidth',3)
Scores2 = [scores; c.scores.ueig]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'y','linewidth',3,'linewidth',3)

% plot([find(Scores == min(Scores)) find(Scores == min(Scores))],[0 max(Scores)],'--','color',[0.5 0.5 0.5])
%plot([length(rlist.edge_indx) length(rlist.edge_indx)],[0 max(Scores)/max(Scores)],'--','color',[0.5 0.5 0.5])
xlim([0 length(Scores2)])
ylim([min(Scores2)/max(Scores2) max(Scores2)/max(Scores2)])
xlabel('Edge')
if strcmp(type,'OD')
    ylabel('OD flow')
elseif strcmp(type,'LargeC')
    ylabel('Largest component')
end
set(gca,'fontsize',16)          % change font size
legend({'closeness','degree','betweenness','eigenvector'},'location','southeast')
legend('boxoff')
set(gca,'fontsize',16)          % change font size
% grid on
h = gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');

print(f,'-dpdf','RECOVERY_TEST_CIR.pdf','-bestfit') % save as pdf file

nrlinks = length(rlist.edge_indx);

temp = c.scores.ucc;
c.R.ucc = nrlinks + 1 - ((temp(1) + temp(end))/2 + sum(temp(2:end-1)))/max(Scores);
temp = c.scores.ud;
c.R.ud = nrlinks + 1 - ((temp(1) + temp(end))/2 + sum(temp(2:end-1)))/max(Scores);
temp = c.scores.ubc;
c.R.ubc = nrlinks + 1 - ((temp(1) + temp(end))/2 + sum(temp(2:end-1)))/max(Scores);
temp = c.scores.ueig;
c.R.ueig = nrlinks + 1 - ((temp(1) + temp(end))/2 + sum(temp(2:end-1)))/max(Scores);


temp = round([ c.R.ucc; c.R.ud; c.R.ubc; c.R.ueig],2);
T = table(temp,'RowNames',{'closeness','degree','betweenness','eigenvector'});






end
