function [KRR_Results] = GSPTB_KernelRidgeRegression_TrainTest(K_train,...
    K_leftout,y_train,y_leftout,Cov_leftout,InnerSplits,PredOpts)

    nTrain = size(K_train,1);
    nTest  = size(K_leftout,1);
    
    assert(size(K_train,2) == nTrain,...
        'Training kernel must be square.');
    
    assert(size(K_leftout,2) == nTrain,...
        'Test kernel must have one column per training subject.');
    
    assert(numel(y_train) == nTrain);
    assert(numel(y_leftout) == nTest);
    assert(numel(Cov_leftout) == nTest);

    for ni = 1:PredOpts.NumInnerFolds
        assert(all(InnerSplits.idx_train{ni} >= 1 & ...
                   InnerSplits.idx_train{ni} <= nTrain));
    
        assert(all(InnerSplits.idx_test{ni} >= 1 & ...
                   InnerSplits.idx_test{ni} <= nTrain));
    
        assert(isempty(intersect( ...
            InnerSplits.idx_train{ni}, ...
            InnerSplits.idx_test{ni})));
    end
    
    y_train   = y_train(:);
    y_leftout = y_leftout(:);
    Cov_leftout = Cov_leftout(:);

    K_outer_train = K_train;

    idx_train = InnerSplits.idx_train;
    idx_test = InnerSplits.idx_test;

    % Initializes the R^2 matrix for the current outer fold
    R2 = NaN(PredOpts.NumInnerFolds,length(PredOpts.Lambda));

    % We now enter the inner loop, where we wish to find the optimal
    % hyperparameter
    for ni = 1:PredOpts.NumInnerFolds

        disp(['Computing inner loop ',num2str(ni),'...']);

        % Training kernel
        K_inner_train = K_outer_train(idx_train{ni},idx_train{ni});
        
        % Test kernel entries; note that we have the similarities
        % between test subjects and training subjects here
        K_inner_test = K_outer_train(idx_test{ni},idx_train{ni});
        
        % Partitioning of behavioural scores
        b_train = y_train(idx_train{ni});
        b_test = y_train(idx_test{ni});

        % Assesses all possible lambda values...
        for l = 1:length(PredOpts.Lambda)

            % The solution is obtained on the training set, for all regularizers
            % Size n_train x n_Lambda, regularized problem by Lambda
            a_train{l} = (K_inner_train + PredOpts.Lambda(l)*eye(size(K_inner_train)))\b_train;

            % Estimated value for the test samples at hand, and the
            % current regularizer
            b_test_fit{ni,l} = K_inner_test*a_train{l};

            R2(ni,l) = GSPTB_ComputeR2(b_test_fit{ni,l},b_test);
        end
    end

    % Finds the parameter value at which we have the best mean R2
    % across inner folds; this is picked as optimum for the current
    % outer fold
    MeanR2 = mean(R2,1,'omitnan');
    [~,idx_best_lambda] = max(MeanR2);
    Lambda_opt = PredOpts.Lambda(idx_best_lambda);

    % Now gets the kernel entries for the left-out data and the whole
    % train/test subset
    K_tt = K_train;
    K_val = K_leftout;

    % Computes the coefficients and subsequent behavioral score
    % estimates
    a_val = (K_tt + Lambda_opt*eye(size(K_tt)))\y_train;
    b_val_fit = K_val*a_val;

    % Evaluates the prediction quality based on R^2, with or without
    % covariate, for the outer fold at hand
    R2_NCV = GSPTB_ComputeR2(b_val_fit,y_leftout);
    R2_NCV_cov = GSPTB_ComputeR2_Covariate(b_val_fit,y_leftout,Cov_leftout);

    disp(['Validation R2 = ',num2str(R2_NCV),'...']);
    disp(['Validation R2 (cova) = ',num2str(R2_NCV_cov),'...']);
    
    KRR_Results.InnerR2 = R2;
    KRR_Results.MeanInnerR2 = MeanR2;
    KRR_Results.y_pred = b_val_fit;
    KRR_Results.R2 = R2_NCV;
    KRR_Results.R2_Cova = R2_NCV_cov;
    KRR_Results.LambdaOpt = Lambda_opt;
end
