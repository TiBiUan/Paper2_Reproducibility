function [InnerSplits] = GSPTB_CreateInnerPredictionSplits(FamilyKernel_Train,PredOpts)

    % We consider the specific set of subjects in the current
    % train/test fold
    G = graph(FamilyKernel_Train);

    % family has size 1 x n_tt_actual, it has one index per family and
    % tags each subject
    family = conncomp(G);

    % Number of families
    families = unique(family);
    F = numel(families);

    % Number of subjects in each family
    FamSize = zeros(F,1);
    
    for f = 1:F
        FamSize(f) = sum(family == families(f));
    end

    % Randomizes within equal/similar sizes, and places larger families first
    % order has size F x 1
    randomTieBreaker = rand(F,1);
    [~,order] = sortrows([-FamSize,randomTieBreaker],[1 2]);

    % We will gradually add families to the smallest current fold, 
    % until all families are assigned. This ensures similarly-sized
    % folds
    FoldSize = zeros(PredOpts.NumInnerFolds,1);
    FamilyFold = zeros(F,1);

    % We go through each family sequentially...
    for fam = 1:F

        % Samples the index of the family at hand
        f = order(fam);
    
        % Puts this family into the currently smallest fold
        [~,k] = min(FoldSize);
        FamilyFold(f) = k;
        FoldSize(k) = FoldSize(k) + FamSize(f);
    end

    % Extracts subject-wise indices linked to the family partitioning
    for k = 1:PredOpts.NumInnerFolds
        isTestFamily = (FamilyFold == k);
        idx_testfams = families(isTestFamily);
        
        idx_test_subspace{k} = find(ismember(family,idx_testfams));
        idx_train_subspace{k} = find(~ismember(family,idx_testfams));

        InnerSplits.idx_test{k} = idx_test_subspace{k};
        InnerSplits.idx_train{k} = idx_train_subspace{k};
    end
end