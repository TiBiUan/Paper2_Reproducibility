function [SparsePredictionResults] = GSPTB_RunSparseActivityPrediction(A_train,A_test,U)
%GSPTB_RUNSPARSEACTIVITYPREDICTION Select a sparsity level on training data and estimate sparse harmonic expression on independent data.
%
%   SparsePredictionResults = GSPTB_RunSparseActivityPrediction(A_train,A_test,U)
%
%   Description
%   -----------
%   Select a sparsity level on training data and estimate sparse harmonic
%   expression on independent data.
%
%   Inputs
%   ------
%   A_train
%       R-by-T1 regional activity matrix used for sparsity-level selection.
%   A_test
%       R-by-T2 independent regional activity matrix used for sparse estimation.
%   U
%       R-by-R harmonic dictionary.
%
%   Outputs
%   -------
%   SparsePredictionResults
%       Structure containing sparse weights, binary activity states, intensity
%       measures, training RMSE curve, test RMSE and R2, and selected sparsity
%       level.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Orthogonal matching pursuit is constrained by a participant- and
%   recording-specific maximum number of active atoms.
%   Positive-only and negative-only activity summaries are returned in nested
%   structures.
%
%   Dependencies
%   ------------
%   GSPTB_FindLOpt
%   mexOMP
%   GSPTB_ComputeMeanIntensity
%   GSPTB_ComputeR2
%
%   See also
%   GSPTB_GFT, GSPTB_ExtractDynamicMetrics


    assert(size(A_train,1) == size(A_test,1), ...
    'Training and testing sessions must contain the same regions.');

    assert(size(U,1) == size(A_train,1), ...
        'The harmonic dictionary and activity data must contain the same regions.');

    % Number of harmonics
    R = size(U,1);
    
    % Creates the dictionary
    D = NaN(R,R);

    for i = 1:size(U,2)
        D(:,i) = U(:,i);
    end

    % Optimal L found on training data
    [L_opt,RMSE_train] = GSPTB_FindLOpt(D,A_train);

    % Now applying to the second set to get the weights on the test run
    param.L = L_opt;
    Weights_test = full(mexOMP(A_test,single(D),param));

    % Selects any non-zero entry as "active"
    BinWeights = (Weights_test ~= 0);

    BinWeights_PosOnly = (Weights_test > 0);
    BinWeights_NegOnly = (Weights_test < 0);

    % These are the main intensity values
    MeanIntensity = GSPTB_ComputeMeanIntensity(BinWeights,Weights_test);

    % These are the values when we only consider positive or negative
    % contributions
    MeanIntensity_PosOnly = GSPTB_ComputeMeanIntensity(BinWeights_PosOnly,Weights_test);
    MeanIntensity_NegOnly = GSPTB_ComputeMeanIntensity(BinWeights_NegOnly,Weights_test);

    % RMSE on the test data for each time point
    A_estimated = double(D)*double(Weights_test);
    RMSE_test = sqrt(sum(((A_test-A_estimated).^2))/R);

    R2_test = NaN(length(RMSE_test),1);

    % R^2 coefficients for each time point
    for t = 1:size(A_test,2)
        R2_test(t) = GSPTB_ComputeR2(A_estimated(:,t),double(A_test(:,t)));
    end

    SparsePredictionResults.Weights = Weights_test;
    SparsePredictionResults.MeanIntensity = MeanIntensity;
    SparsePredictionResults.BinWeights = BinWeights;
    SparsePredictionResults.RMSE_train = RMSE_train';
    SparsePredictionResults.RMSE = RMSE_test';
    SparsePredictionResults.R2 = R2_test;
    SparsePredictionResults.L_opt = L_opt;

    SparsePredictionResults.PosOnly.MeanIntensity = MeanIntensity_PosOnly;
    SparsePredictionResults.PosOnly.BinWeights = BinWeights_PosOnly;

    SparsePredictionResults.NegOnly.MeanIntensity = MeanIntensity_NegOnly;
    SparsePredictionResults.NegOnly.BinWeights = BinWeights_NegOnly;
end
