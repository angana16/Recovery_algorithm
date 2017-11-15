% Input: 
% G:        original graph
% OD:       Original OD matrix
% rlist:    list of removed edges
% budget:   number of edges that can be restored
% Output:
% sset:     order of edge recovery
% scores:   scores
% evalNume: number of evaluations in each iteration
function [sset,scores,evalNum] = greedy_lazy(G, OD, rlist, budget,type)

n = size(G.Edges,1);                    % number of edges
deltas = inf*ones(1,length(rlist.edge_indx)); % initialize 
eps = 1e-3;                             % tolerance
sset = [];                              % start with an empty or specified set
% get the OD flow of the disrupted network
Gr = G;                 % original graph
Gr = rmedge(Gr,rlist.edge_indx);  % remove the edges in the rlist
test = [];
curVal = inc_recover(Gr,[],[],test,OD,type); % current value
% get current OD flow
evalNum = zeros(1,n); scores = [];      % keep track of statistics

i = 0;
while 1
    if i > budget
        break;
    end
    i = i+1;
    bestimprov = 0;
    evalNum(i) = 0;
    
    [~,order] = sort(deltas,'descend');
    
    for test = order
            evalNum(i) = evalNum(i) + 1;
            % get improvement in solution
            improv = inc_recover(Gr,sset,test,rlist.nodes_indx(test,1:2),OD,type) - curVal;
            deltas(test) = improv;
            bestimprov = max(bestimprov,improv);
    end
    % pick the solution that gives the best improvement
    %-----------------------
    % could be more than one
    % pick the first one
%     argmax = find(deltas==max(deltas),1); 
    % randomly pick
    argmax = find(deltas==max(deltas));
    argmax = argmax(randi(length(argmax),1));
        
    if deltas(argmax) > eps % nontrivial improvement by adding argmax
        sset = [sset,argmax];       % sset is the index of the edge in rlist NOT the index of the edge in G
        curVal = curVal + deltas(argmax);
        scores(i) = curVal;
        % update G
        Gr = addedge(Gr,rlist.nodes_indx(argmax,1),rlist.nodes_indx(argmax,2),1);
    else
        break
    end
end
