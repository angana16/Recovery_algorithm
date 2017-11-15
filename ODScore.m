function newScore = ODScore(G,OD,type)
% get the OD score of the network

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