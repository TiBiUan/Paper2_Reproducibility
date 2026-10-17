function [KRR_Results] = GSPTB_KernelRidgeRegression(Ker,b,PredOpts)
%GSPTB_KERNELRIDGEREGRESSION Perform nested cross-validated kernel ridge regression.
%
%   KRR_Results = GSPTB_KernelRidgeRegression(Ker,b,PredOpts)
%
%   Description
%   -----------
%   Perform nested cross-validated kernel ridge regression.
%
%   Inputs
%   ------
%   Ker
%       S-by-S subject kernel matrix.
%   b
%       S-element target vector.
%   PredOpts
%       Prediction-options structure defining splits, family constraints,
%       candidate regularization values, covariates, and nested cross-validation
%       dimensions.
%
%   Outputs
%   -------
%   KRR_Results
%       Structure containing outer-fold R2 values, covariate-adjusted R2 values,
%       and selected regularization parameters.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Family-related subjects are kept together when random splits are
%   generated.
%   If all split fields are populated, predefined outer and inner partitions
%   are used.
%   Random splitting depends on the active MATLAB random-number generator
%   state.
%
%   Dependencies
%   ------------
%   GSPTB_CreateSplit
%   GSPTB_ComputeR2
%   GSPTB_ComputeR2_Covariate
%
%   See also
%   GSPTB_ComputeHarmonicRangeKernel


    idx_leftout = PredOpts.Splits.leftout;
    idx_tt = PredOpts.Splits.tt;
    idx_test = PredOpts.Splits.test;
    idx_train = PredOpts.Splits.train;

    split_inputs = {idx_leftout,idx_tt,idx_test,idx_train};
    is_empty = cellfun(@isempty,split_inputs);
    
    if all(is_empty)
        is_predef = false;
    elseif ~any(is_empty)
        is_predef = true;
    else
        error('Predefined splitting is incomplete!');
    end

    if is_predef
        assert(numel(idx_tt) == PredOpts.NumOuterSplits, ...
            'idx_tt must contain n_outer entries!');
    
        assert(numel(idx_leftout) == PredOpts.NumOuterSplits, ...
            'idx_leftout must contain n_outer entries!');
    
        for no = 1:PredOpts.NumOuterSplits
            assert(numel(idx_train{no}) == PredOpts.NumInnerFolds, ...
                'idx_train{%d} must contain KCV folds!',no);
    
            assert(numel(idx_test{no}) == PredOpts.NumInnerFolds, ...
                'idx_test{%d} must contain KCV folds!',no);
        end
    end

    % Number of subjects
    S = size(Ker,2);

    assert(size(Ker,1)==S,'Kernel matrix must be square!');
            
    % Left out subjects (outer loop) used, at the end, to ascertain the
    % prediction quality
    n_leftout = floor(PredOpts.LeftOutFraction*S);

    % Subjects for training and assessing hyperparameter
    n_tt = S - n_leftout;

    % Outer loop
    for no = 1:PredOpts.NumOuterSplits

        disp(['Computing outer loop ',num2str(no),'...']);

        % If we have not provided splitting information, we compute it
        if ~is_predef

            % We randomly sample subjects to include in either set, while
            % acounting for families
            % Note that the exact numbers inside idx_tt and idx_leftout may 
            % not exactly match the specified ones, due to family constraints
            [idx_tt{no},idx_leftout{no}] = GSPTB_CreateSplit(n_tt,n_leftout,PredOpts.FamilyKernel);
    
            % Samples the family relationships inside the train/test only
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
                
                idx_test_subspace{no}{k} = find(ismember(family,idx_testfams));
                idx_train_subspace{no}{k} = find(~ismember(family,idx_testfams));

                idx_test{no}{k} = idx_tt{no}(idx_test_subspace{no}{k});
                idx_train{no}{k} = idx_tt{no}(idx_train_subspace{no}{k});
            end
        end

        % Our behavioural scores are also extracted similarly, so are
        % covariate values
        b_tt{no} = b(idx_tt{no});
        b_val{no} = b(idx_leftout{no});
        cova_val{no} = PredOpts.Covariate(idx_leftout{no},:);

        % Initializes the R^2 matrix for the current outer fold
        R2 = NaN(PredOpts.NumInnerFolds,length(PredOpts.Lambda));

        % We now enter the inner loop, where we wish to find the optimal
        % hyperparameter
        for ni = 1:PredOpts.NumInnerFolds

            disp(['Computing inner loop ',num2str(ni),'...']);

            % Training kernel
            K_train = Ker(idx_train{no}{ni},idx_train{no}{ni});
            
            % Test kernel entries; note that we have the similarities
            % between test subjects and training subjects here
            K_test = Ker(idx_test{no}{ni},idx_train{no}{ni});
            
            % Partitioning of behavioural scores
            b_train = b(idx_train{no}{ni});
            b_test = b(idx_test{no}{ni});

            % Assesses all possible lambda values...
            for l = 1:length(PredOpts.Lambda)

                % The solution is obtained on the training set, for all regularizers
                % Size n_train x n_Lambda, regularized problem by Lambda
                a_train{l} = (K_train + PredOpts.Lambda(l)*eye(size(K_train)))\b_train;

                % Estimated value for the test samples at hand, and the
                % current regularizer
                b_test_fit{ni,l} = K_test*a_train{l};

                R2(ni,l) = GSPTB_ComputeR2(b_test_fit{ni,l},b_test);
            end
        end

        % Finds the parameter value at which we have the best mean R2
        % across inner folds; this is picked as optimum for the current
        % outer fold
        tmp = find(mean(R2) == max(mean(R2)));
        idx_best_lambda(no) = tmp(1);
        Lambda_opt(no) = PredOpts.Lambda(idx_best_lambda(no));

        % Now gets the kernel entries for the left-out data and the whole
        % train/test subset
        K_tt = Ker(idx_tt{no},idx_tt{no});
        K_val = Ker(idx_leftout{no},idx_tt{no});
    
        % Computes the coefficients and subsequent behavioral score
        % estimates
        a_val{no} = (K_tt + Lambda_opt(no)*eye(size(K_tt)))\b_tt{no};
        b_val_fit{no} = K_val*a_val{no};

        % Evaluates the prediction quality based on R^2, with or without
        % covariate, for the outer fold at hand
        R2_NCV(no) = GSPTB_ComputeR2(b_val_fit{no},b_val{no});
        R2_NCV_cov(no) = GSPTB_ComputeR2_Covariate(b_val_fit{no},b_val{no},cova_val{no});

        disp(['Validation R2 for fold ',num2str(no),' = ',num2str(R2_NCV(no)),'...']);
        disp(['Validation R2 (cova) for fold ',num2str(no),' = ',num2str(R2_NCV_cov(no)),'...']);
    end

    KRR_Results.R2 = R2_NCV;
    KRR_Results.R2_Cova = R2_NCV_cov;
    KRR_Results.LambdaOpt = Lambda_opt;
end
