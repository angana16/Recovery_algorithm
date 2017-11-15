clc
clear
close all hidden

% build mbta network

% load data
%------------------------------------------------
load('mbta_data.mat')
% load('mbta_data1.mat')  % with winter 2015 event
% includes adjacency matrix, list of nodes, and list of nodes impacted by
% snow storm

G = graph(Adj);
G.Nodes.id = nodes.id;
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
temp = 5*ones(nnodes,1);
temp(nodes_snow.index) = 5;
G.Nodes.Snow = temp;
p.MarkerSize = temp;
p.XData = nodes.lat;
p.YData = nodes.lon;

% labelnode(p,nodes_snow.index,nodes_snow.id)
set(gca,'fontsize',16)          % change font size
% grid on
f = gcf;
set(f,'PaperPositionMode','auto');         
set(f,'PaperOrientation','landscape');
set(gca, 'XTick', [],'YTick', [],'Box', 'off','XColor','none','YColor','none','Color','none');
print(f,'-dpdf','fig_mbta_layout.pdf','-bestfit') % save as pdf file
print(f,'-dpng','fig_mbta_layout.png') % save as pdf file

% G.Nodes.NodeColors = degree(G);     % set nodes colors equal to node degree
% p.NodeCData = G.Nodes.NodeColors;   % set
% p.LineWidth = randi([1 5],size(G.Edges,1),1);

% load and plot the power network
draw = 'false';
if strcmp(draw,'true')
    load('power_adj.mat')
    Gp = graph(adjpower);
    nnodes = size(Gp.Nodes,1);
    nlinks = size(Gp.Edges,1);
    f = figure;
    p = plot(Gp,'layout','force');
    p.Marker = 'o';
    p.NodeColor = [0.5430 0 0] ;
    p.EdgeColor = [0.5430 0 0] ;
    temp = 5*ones(nnodes,1);
    p.MarkerSize = temp;
    set(gca,'fontsize',16)          % change font size
    % grid on
    f = gcf;
    set(f,'PaperPositionMode','auto');
    set(f,'PaperOrientation','landscape');
    set(gca, 'XTick', [],'YTick', [],'Box', 'off','XColor','none','YColor','none','Color','none');
    print(f,'-dpdf','fig_power_layout.pdf','-bestfit') % save as pdf file
    print(f,'-dpng','fig_power_layout.png') % save as pdf file
    
end


% made up OD matrix
%------------------------------------------------
OD = 10*(ones(nnodes) - diag(ones(nnodes,1)));

% list of edges that are initially removed from the network
%------------------------------------------------
% rlist.edge_indx = [1:3 30:35];            % index of edges that are removed from the full graph
var1 = nodes_snow.index;                    % nodes based on input data
list.nodes_indx = table2array(G.Edges);     % get edge start/end nodes
list.nodes_indx(:,end) = [];                % remove weights
var2 = ismember(list.nodes_indx,var1);
var3 = +(sum(var2,2)>0);
rlist.edge_indx = find(var3);
rlist.nodes_indx = table2array(G.Edges(rlist.edge_indx,:));
rlist.nodes_indx(:,end) = [];
budget = length(rlist.edge_indx);        % available budget

% plot removed edges
%------------------------------------------------
Gr = graph(Adj);
Gr = rmedge(Gr,rlist.edge_indx);
draw = 'true';
if strcmp(draw,'true')
    f = figure;
    subplot(1,2,1)
    p = plot(G,'layout','force');
    p.XData = nodes.lat;
    p.YData = nodes.lon;

    subplot(1,2,2)
    pr = plot(Gr,'layout','force');
    pr.XData = p.XData;
    pr.YData = p.YData;
end
%------------------------------------------------


%------------------------------------------------
% pick the type of objective function to optimize
% type = 'OD';        % od matrix
type = 'LargeC';  % largest component
%------------------------------------------------

% greedy recovery
%------------------------------------------------
% G - original graph; OD - matrix; rlist - list of removed edgest; budget;
% type of functionality
[greedy.sset,greedy.scores,greedy.evalNum] = greedy_lazy(G, OD, rlist, budget,type);

%------------------------------------------------

% plot failure & recovery curve
%------------------------------------------------
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

draw = 'false';
if strcmp(draw,'true')
    f = figure();
    plot(Scores/max(Scores))
    hold on
    % plot([find(Scores == min(Scores)) find(Scores == min(Scores))],[0 max(Scores)],'--','color',[0.5 0.5 0.5])
    plot([length(rlist.edge_indx) length(rlist.edge_indx)],[0 max(Scores)/max(Scores)],'--','color',[0.5 0.5 0.5])
    xlim([0 length(Scores)])
    ylim([0 max(Scores)/max(Scores)])
    xlabel('Edge')
    if strcmp(type,'OD')
        ylabel('OD flow')
    elseif strcmp(type,'LargeC')
        ylabel('Largest component')
    end
    set(gca,'fontsize',16)          % change font size
end

% nodal centrality
%------------------------------------------------
c.ucc = centrality(G,'closeness');
c.ud = centrality(G,'degree');
c.ubc = centrality(G,'betweenness');
c.ueig = centrality(G,'eigenvector');

draw = 'false';
if strcmp(draw,'true')
    % plot distribution of centrality measures
    %------------------------------------------------
    f = figure;
    p = plot(G,'layout','force');
    p.XData = nodes.lat;
    p.YData = nodes.lon;
    p.NodeCData = c.ucc;
    colormap jet
    colorbar
    title('Closeness Centrality Scores - Unweighted')
    
    f = figure;
    p = plot(G,'layout','force');
    p.XData = nodes.lat;
    p.YData = nodes.lon;
    p.NodeCData = c.ubc;
    colormap jet
    colorbar
    title('Betweenness Centrality Scores - Unweighted')
end


% recovery: greedy, ucc, ud, ubc, ueig
%------------------------------------------------

rnodes_indx = nodes_snow.index;        % removed nodes index
% closenesspip in
[c.scores.ucc,c.reclist.ucc] = recover_centrality(G,OD,c.ucc,type,rnodes_indx,list);
[c.scores.ud,c.reclist.ud] = recover_centrality(G,OD,c.ud,type,rnodes_indx,list);
[c.scores.ubc,c.reclist.ubc] = recover_centrality(G,OD,c.ubc,type,rnodes_indx,list);
[c.scores.ueig,c.reclist.ueig] = recover_centrality(G,OD,c.ueig,type,rnodes_indx,list);


%------------------------------------------------
% call cross-entropy algorithm
CE_v1
%------------------------------------------------


% plot
%------------------------------------------------
f = figure();
plot(Scores/max(Scores),'b','linewidth',3)   % greedy
hold on
Scores2 = [scores; c.scores.ucc]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'m','linewidth',3)
Scores2 = [scores; c.scores.ud]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'k','linewidth',3)
Scores2 = [scores; c.scores.ubc]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'g','linewidth',3)
Scores2 = [scores; c.scores.ueig]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'y','linewidth',3,'linewidth',3)

Scores2 = [scores; ce.scores]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'r','linewidth',3)

plot(scores/max(Scores2),'color',[0.5 0.5 0.5],'linewidth',3)
% plot([find(Scores == min(Scores)) find(Scores == min(Scores))],[0 max(Scores)],'--','color',[0.5 0.5 0.5])
plot([length(rlist.edge_indx) length(rlist.edge_indx)],[0 max(Scores)/max(Scores)],'--','color',[0.5 0.5 0.5])
xlim([0 length(Scores)])
ylim([0 max(Scores)/max(Scores)])
xlabel('Edge')
if strcmp(type,'OD')
    ylabel('OD flow')
elseif strcmp(type,'LargeC')
    ylabel('Largest component')
end
set(gca,'fontsize',16)          % change font size
legend({'greedy','closeness','degree','betweenness','eigenvector','cross entropy'},'location','southeast')
legend('boxoff')
set(gca,'fontsize',16)          % change font size
% grid on
h = gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
print(f,'-dpdf','fig_mbta_scores.pdf','-bestfit') % save as pdf file


% get final score: area under the nominal value
nrlinks = length(rlist.edge_indx);
temp = greedy.scores;
greedy.R = nrlinks + 1 - ((temp(1) + temp(end))/2 + sum(temp(2:end-1)))/max(Scores);
temp = c.scores.ucc;
c.R.ucc = nrlinks + 1 - ((temp(1) + temp(end))/2 + sum(temp(2:end-1)))/max(Scores);
temp = c.scores.ud;
c.R.ud = nrlinks + 1 - ((temp(1) + temp(end))/2 + sum(temp(2:end-1)))/max(Scores);
temp = c.scores.ubc;
c.R.ubc = nrlinks + 1 - ((temp(1) + temp(end))/2 + sum(temp(2:end-1)))/max(Scores);
temp = c.scores.ueig;
c.R.ueig = nrlinks + 1 - ((temp(1) + temp(end))/2 + sum(temp(2:end-1)))/max(Scores);
temp = ce.scores;
ce.R = nrlinks + 1 - ((temp(1) + temp(end))/2 + sum(temp(2:end-1)))/max(Scores);


temp = round([greedy.R; c.R.ucc; c.R.ud; c.R.ubc; c.R.ueig; ce.R],2);
T = table(temp,'RowNames',{'greedy','closeness','degree','betweenness','eigenvector','cross entropy'})

