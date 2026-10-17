function [PredSplits] = GSPTB_CreatePredictionSplits(PredOpts)

    % Number of subjects
    S = size(PredOpts.FamilyKernel,1);
    
    % How many subjects in each side of the outer split
    n_leftout = floor(PredOpts.LeftOutFraction*S);
    n_tt = S - n_leftout;

    % Splitting for each outer fold
    for no = 1:PredOpts.NumOuterSplits
    
        % Splitting (outer loop)
        [idx_tt{no},idx_leftout{no}] = GSPTB_CreateSplit(n_tt,n_leftout,PredOpts.FamilyKernel);
        
        % Samples the Family relationships
        K_Family_tt = PredOpts.FamilyKernel(idx_tt{no},idx_tt{no});
        
        % We consider the specific set of subjects in the current
        % train/test fold
        G = graph(K_Family_tt);
        
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
            
            idx_test{no}{k} = find(ismember(family,idx_testfams));
            idx_test_fullsubj{no}{k} = idx_tt{no}(idx_test{no}{k});

            idx_train{no}{k} = find(~ismember(family,idx_testfams));
            idx_train_fullsubj{no}{k} = idx_tt{no}(idx_train{no}{k});
        end  
    end

    % Summarizes the splitting information to return
    PredSplits.leftout = idx_leftout;
    PredSplits.tt = idx_tt;
    PredSplits.test = idx_test_fullsubj;
    PredSplits.train = idx_train_fullsubj;
end