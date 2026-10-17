function [idx1,idx2,le,le2] = GSPTB_CreateSplit(n1,n2,K_Family)

% This contains all our indices, which we want to partition
idx = randperm(n1+n2);
idx_fin = idx;

% Determines complete family relationships
G = graph(K_Family);
family = conncomp(G);

% This will contain the first partition
idx1 = [];

% We gradually append complete families until the desired size is reached
i = 1;

while length(idx1) < n1

    % Samples the next candidate subject
    tmp_idx = idx(1);

    % Finds all subjects belonging to the same family component
    tmp_fam = family(tmp_idx);
    tmp = find(family == tmp_fam);

    % Adds the complete family
    idx1 = [idx1,tmp];
    idx1 = unique(idx1);

    le(i) = length(idx1);

    % Removes the complete family from the remaining candidates
    idx(ismember(idx,tmp)) = [];

    le2(i) = length(idx);

    i = i+1;
end

% Remaining subjects form the second partition
idx2 = idx_fin(~ismember(idx_fin,idx1));

end