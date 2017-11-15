clc
clear
close all hidden

% build mbta network

% load data
%------------------------------------------------
load('Indian_datafiles.mat')
% load('mbta_data1.mat')  % with winter 2015 event
% includes adjacency matrix, list of nodes, and list of nodes impacted by
% snow storm
Adj = adj_IRN;
nodes = nodes_irn;
G = graph(adj_IRN);
G.Nodes.id = nodes_irn.id;
nnodes = size(G.Nodes,1);
nlinks = size(G.Edges,1);

f = figure;
%worldmap('india')

p = plot(G,'layout','force');
p.Marker = 'o';
p.NodeColor = [0 0 0];
% p.NodeFaceColor = [1 1 1];
p.EdgeColor = [0 0 0];
temp = 5*ones(nnodes,1);
temp(tsunami.index) = 8;
G.Nodes.tsunami = temp;
p.MarkerSize = temp;


temp = zeros(nnodes,3);
for i=1:length(tsunami.index)
    temp(tsunami.index(i),:) = [1 0 0];
end 
p.NodeColor = temp;
p.YData = nodes_irn.lat;
p.XData = nodes_irn.lon;

% labelnode(p,nodes_snow.index,nodes_snow.id)
set(gca,'fontsize',16)          % change font size
% grid on
f = gcf;
set(f,'PaperPositionMode','auto');         
set(f,'PaperOrientation','landscape');
set(gca, 'XTick', [],'YTick', [],'Box', 'off','XColor','none','YColor','none','Color','none');
print(f,'-dpdf','fig_india_layout.pdf','-bestfit') % save as pdf file
print(f,'-dpng','fig_india_layout.png') % save as pdf file


var1 = tsunami.index;                    % nodes based on input data
list.nodes_indx = table2array(G.Edges);     % get edge start/end nodes
list.nodes_indx(:,end) = [];                % remove weights
var2 = ismember(list.nodes_indx,var1);
var3 = +(sum(var2,2)>0);
rlist.edge_indx = find(var3);
rlist.nodes_indx = table2array(G.Edges(rlist.edge_indx,:));
rlist.nodes_indx(:,end) = [];
budget = length(rlist.edge_indx);

%------------------------------------------------

% made up OD matrix
%------------------------------------------------
OD = adj_IRN_weighted;

% list of edges that are removed from tsunami
%------------------------------------------------
% rlist.edge_indx = [1:3 30:35];            % index of edges that are removed from the full graph
var1 = tsunami.index;                    % nodes based on input data
list.nodes_indx = table2array(G.Edges);     % get edge start/end nodes
list.nodes_indx(:,end) = [];                % remove weights
var2 = ismember(list.nodes_indx,var1);
var3 = +(sum(var2,2)>0);
rlist.edge_indx = find(var3);
rlist.nodes_indx = table2array(G.Edges(rlist.edge_indx,:));
rlist.nodes_indx(:,end) = [];
budget = length(rlist.edge_indx);        % available budget



%------------------------------------------------
% pick the type of objective function to optimize
 type = 'OD';        % od matrix
%type = 'LargeC';  % largest component
%------------------------------------------------
% greedy recovery
%------------------------------------------------
% G - original graph; OD - matrix; rlist - list of removed edgest; budget;
% type of functionality
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

Scores = [scores; greedy.scores']; % combine failure and recovery scores
%% 

draw = 'true';
if strcmp(draw,'true')
    f = figure();
    plot(Scores/max(Scores))
    hold on
    % plot([find(Scores == min(Scores)) find(Scores == min(Scores))],[0 max(Scores)],'--','color',[0.5 0.5 0.5])
    plot([length(rlist.edge_indx) length(rlist.edge_indx)],[0 max(Scores)/max(Scores)],'--','color',[0.5 0.5 0.5])
    xlim([0 length(Scores)])
    ylim([0.8 max(Scores)/max(Scores)])
    xlabel('Edge')
    if strcmp(type,'OD')
        ylabel('OD flow')
    elseif strcmp(type,'LargeC')
        ylabel('Largest component')
    end
    set(gca,'fontsize',16)          % change font size
end

c.ucc = centrality(G,'closeness');
c.ud = centrality(G,'degree');
c.ubc = centrality(G,'betweenness');
c.ubc = 2*c.ubc/((nnodes-2)*(nnodes-1));
maxim = find(c.ubc == max(c.ubc(:)));
c.ueig = centrality(G,'eigenvector');


draw = 'false';
if strcmp(draw,'true')
    % plot distribution of centrality measures
    %------------------------------------------------
    f = figure;
    p = plot(G,'layout','force');
    p.XData = nodes.lon;
    p.YData = nodes.lat;
    p.NodeCData = c.ucc;
    colormap jet
    colorbar
    title('Closeness Centrality Scores - Unweighted')
    
    f = figure;
    p = plot(G,'layout','force');
    p.XData = nodes.lon;
    p.YData = nodes.lat;
    p.NodeCData = c.ubc;
    colormap jet
    colorbar
    title('Betweenness Centrality Scores - Unweighted')
end

% recovery: greedy, ucc, ud, ubc, ueig
%------------------------------------------------

rnodes_indx = tsunami.index;        % removed nodes index
% closeness
[c.scores.ucc,c.reclist.ucc] = recover_centrality(G,OD,c.ucc,type,rnodes_indx,list);
[c.scores.ud,c.reclist.ud] = recover_centrality(G,OD,c.ud,type,rnodes_indx,list);
[c.scores.ubc,c.reclist.ubc] = recover_centrality(G,OD,c.ubc,type,rnodes_indx,list);
[c.scores.ueig,c.reclist.ueig] = recover_centrality(G,OD,c.ueig,type,rnodes_indx,list);

CE_v1
