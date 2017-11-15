% mbta visualizations

clc
clear 
close all hidden

% load data
%------------------------------------------------
load('mbta_data.mat')       % mbta data
load('mbta_ce3.mat')        % simulation results
load('edges_recovery.mat')  % recovery edges for plotting purposes

% plot
%------------------------------------------------
f = figure();
plot(Scores/max(Scores),'b','linewidth',2)   % greedy
hold on
Scores2 = [scores; c.scores.ucc]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'m','linewidth',2)
Scores2 = [scores; c.scores.ud]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'k','linewidth',2)
Scores2 = [scores; c.scores.ubc]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'g','linewidth',2)
Scores2 = [scores; c.scores.ueig]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'y','linewidth',2)

Scores2 = [scores; ce.scores]; % combine failure and recovery scores
plot(Scores2/max(Scores2),'r','linewidth',2)

plot(scores/max(Scores2),'color',[0.5 0.5 0.5],'linewidth',2)
% plot([find(Scores == min(Scores)) find(Scores == min(Scores))],[0 max(Scores)],'--','color',[0.5 0.5 0.5])
plot([length(rlist.edge_indx) length(rlist.edge_indx)],[0 max(Scores)/max(Scores)],'--','color',[0.5 0.5 0.5])
xlim([0 length(Scores)])
ylim([0 max(Scores)/max(Scores)])
xlabel('System state')
if strcmp(type,'OD')
    ylabel('OD flow')
elseif strcmp(type,'LargeC')
    ylabel('System functionality')
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

%------------------------------------------
% graphs
%------------------------------------------
% gif file to illustrate network recovery

%------------------------------------------
Gr = graph(Adj);
Gr = rmedge(Gr,rlist.edge_indx);
scnsize = get(0,'ScreenSize');
delayt = 1.5;

% select output file
%--------------------------
% outfile = 'gif_mbta_recovery_greedy.gif';     % greedy
outfile = 'gif_mbta_recovery_ce.gif';           % ce
% outfile = 'gif_mbta_recovery_ubc.gif';        % ubc


% get the list and order of recovery edges and their start & end nodes
% select solution
%--------------------------
var1 = table2array(G.Edges(ce.edges,:));                    % ce
% var1 = table2array(G.Edges(edges_recovery.ubc,:));        % ubc
% var1 = table2array(G.Edges(edges_recovery.greedy,:));     % greedy
var1(:,end) = [];

for k = 1:size(var1,1)
    close all
    f = figure('Position',scnsize);
    % increment the graph by adding edges
    Gr = addedge(Gr,var1(k,1),var1(k,2),2);
    p = plot(Gr,'layout','force');
    p.MarkerSize = 6;
    p.LineWidth = 3;
    p.Marker = 'o';
    p.NodeColor = [0 0 0];
    temp = find(Gr.Edges.Weight==2);
    tempcolor = zeros(size(Gr.Edges,1),3);
    tempcolor(temp,1) = 1;
    %         p.EdgeColor = [0 0 0];
    tempwidth = 2*ones(size(Gr.Edges,1),1);
    tempwidth(temp) = 7;
    p.LineWidth = tempwidth;
    p.EdgeColor = tempcolor;
    %         set(gca, 'XTick', [],'YTick', [],'Box', 'off','XColor','none','YColor','none','Color','none');
    set(gcf,'Color',[1 1 1])
    set(gca,'xtick',[],'ytick',[])
    box off
    axis off
    axis equal
    p.XData = nodes.lat;
    p.YData = nodes.lon;
    xlim([0.999*min(nodes.lat) 1.001*max(nodes.lat)])
    ylim([min(nodes.lon) max(nodes.lon)])
    
    % gif utilities
    set(gcf,'color','w');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % On first iteration, create file. Otherwise, append to file.
    if k == 1
        %         imwrite(imind,cm,outfile,'gif','DelayTime',delayt,'loopcount',0);
        %         imwrite(imind,cm,outfile,'gif','DelayTime',delayt,'loopcount',inf);
        imwrite(imind,cm,outfile,'gif','DelayTime',delayt,'loopcount',1);
    else
        imwrite(imind,cm,outfile,'gif','DelayTime',delayt,'writemode','append');
    end
end



%------------------------------------------
% scores
%------------------------------------------
% gif file to illustrate functionality recovery

scnsize = get(0,'ScreenSize');
delayt = 1.5;
outfile = 'gif_mbta_recovery_scores.gif';

% get the list and order of recovery edges and their start & end nodes
var1 = ce.scores;       % ce
var2 = c.scores.ubc;    % ubc

for k = 1:size(var1,1)
    close all
%     f = figure('Position',scnsize);
        f = figure();

    plot(scores/max(Scores2),'color',[0.5 0.5 0.5],'linewidth',7)

    hold on
    temp = c.scores.ubc(1:k);
    plot(35 : 35 + k, [scores(end)/max(Scores2) temp'/max(Scores2)],'-g','linewidth',7)
    plot(35 + k, temp(end)/max(Scores2),'go','markerfacecolor','g','markersize',20)
    
    temp = ce.scores(1:k);
    plot(35 : 35 + k, [scores(end)/max(Scores2) temp'/max(Scores2)],'-r','linewidth',7)
    plot(35 + k, temp(end)/max(Scores2),'rs','markerfacecolor','r','markersize',20)

    xlim([0 length(Scores)])
    ylim([0 max(Scores)/max(Scores)])
    xlabel('System state')
    ylabel('System functionality')
    set(gca,'fontsize',22)          % change font size
    set(gcf,'Color',[1 1 1])
%     set(gca,'xtick',[])
%     box off
%     axis off
    %         axis equal
    % gif utilities
    set(gcf,'color','w');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % On first iteration, create file. Otherwise, append to file.
    if k == 1
        imwrite(imind,cm,outfile,'gif','DelayTime',delayt,'loopcount',1);
    else
        imwrite(imind,cm,outfile,'gif','DelayTime',delayt,'writemode','append');
    end
end


