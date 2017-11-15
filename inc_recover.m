function newScore = inc_recover(G,sset,test,testnodes,OD,type)

% find start and end nodes of an edge

if ~isempty(test)
    if ~ismember(test,sset)
        G = addedge(G,testnodes(1),testnodes(2),1);
    end
end


bins = conncomp(G);
temp = unique(bins);
flow = 0;
sizeC = 0;
for i = 1:length(temp)
    indx = find(bins == temp(i));
    ODsub = OD(indx,indx);
    flow = sum(sum(triu(ODsub))) + flow;
    sizeC = max(sizeC,length(indx));
end
if strcmp(type,'OD')
    newScore = flow;
elseif strcmp(type,'LargeC')
    newScore = sizeC;
end
