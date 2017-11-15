function [S] = sample_discrete(p0,nsample)

S = zeros(nsample,size(p0,2));
parfor j = 1:nsample
    sample = [];
    if find(sum(p0,2) > 1+1e-3)
        error('not valid probability')
    end
    x = discretesample(p0(1,:), 1); % first sample
    sample(1) = x;
    for i = 1:size(p0,1)-1    % no repetition % once an edge is sampled it cannot be sampled again
        % update sampling
        var1 = p0(i+1,:);
        var2 = sum(var1(sample));
        var1([1:x-1 x+1:end]) = var1([1:x-1 x+1:end]) + var2/(length(var1) - length(sample));
        var1(sample) = 0;
        x = discretesample(var1, 1);
        sample = [sample x];
    end
    S(j,:) = sample;
end
end
