%% 1. Behavioral prediction
%
% Cognition is predicted using kernel ridge regression within a 10-fold
% family-aware cross-validation scheme.
%
% For each outer fold:
%
%   1. Population connectome harmonics are estimated from training
%      subjects only.
%   2. Graph-frequency cutoffs are estimated from training subjects only.
%   3. Prediction features are computed using the resulting population
%      representation and subject-specific representations.
%   4. Train and test kernels are constructed separately.
%   5. Hyperparameters are selected using family-aware inner folds
%      restricted to the outer training population.
%   6. Predictions are generated for the held-out subjects.
%
% Sixteen feature sets are evaluated:
%
% Unfiltered:
%   FC                  Functional connectivity
%   Cx                  Second-order moment
%
% Filter-dependent:
%   DC                  Difference in second-order moments
%   EC                  Error second-order moment
%   CT                  Cross-term
%   FC_Pop              FC after population-based filtering
%   FC_Ind              FC after subject-specific filtering
%   Delta_FC            Difference in FC
%   Error_FC            FC of the functional error
%
% Filter-dependent quantities are evaluated separately for Low and High
% graph frequencies
%
% Computational time: around 20 minutes
disp('Starting to run Practical Utility analyses (part 3)...');

disp('1. Behavioral prediction...');


% Configuration
rng(1);

% Family-aware outer prediction splits
PredOptsSplit.FamilyKernel = Data.FamilyKernel;
PredOptsSplit.NumFolds = 10;

% Kernel ridge regression parameters
PredOpts.NumInnerFolds = Cfg.Prediction.NumInnerFolds;
PredOpts.Lambda = Cfg.Prediction.Lambda;

% Creates outer prediction folds
% This splitting differs from the one used elsewhere in the manuscript, as
% we want to ensure each subject is present in a left-out set only once
PredSplits = GSPTB_CreateBehavioralPredictionSplits(PredOptsSplit);

% Functional data parameters
id_REST1 = 8;
id_session = 1;

% Filter sub-types
Filters = {'Low','High'};

% Handy parameters

% Number of subjects
S = numel(Data.Adj);

% Number of harmonics
nHarmonics = size(Data.Adj{1},1);

% Number of edges without the diagonal
nEdges = nHarmonics*(nHarmonics-1)/2;

% Number of edges with the diagonal
nEdgesUpper = nHarmonics*(nHarmonics+1)/2;

% Loads resting-state data
X_OI = GSPTB_LoadTaskData(fullfile(DataDir,'Inputs'),'REST1',...
    Cfg.Atlas.CorticalIndices,Cfg.Atlas.IsCortexOnly);

% Makes sure we do not use a previously defined G (for groupings)
clear G

% Names of all the quantities to assess predictive potential for, including
% Functional connectivity and second-order moment (no filtering)
% Difference in second-order moment matrices for either filtering
% Decomposition into error and cross-term for either filtering
% Population-filtering (de)coupled FC
% Individual-filtering (de)coupled FC
% Difference in functional connectivity (pop - ind)
% Functional connectivity of the error signal for either filtering
FeatureNames = {'Cx','FC',...
    'DC_Low','DC_High',...
    'EC_Low','EC_High','CT_Low','CT_High',...
    'FC_Pop_Low','FC_Pop_High',...
    'FC_Ind_Low','FC_Ind_High',...
    'Delta_FC_Low','Delta_FC_High',...
    'Error_FC_Low','Error_FC_High'};

% Runs the same for each outer prediction folds
for k = 1:PredSplits.NumFolds

    k

    % Training and test populations
    idPop = PredSplits.Fold(k).Train;
    idTest = PredSplits.Fold(k).Test;

    % Training-derived population representations
    A_Pop = GSPTB_GetPopulationAdjacencyMatrices(Data,Cfg,idPop);

    % For the AllEdges/AllValues case, we get
    % harmonics/eigenvalues/Laplacian
    H_Pop = GSPTB_GetPopulationConnectomeHarmonics(A_Pop.AllEdges.AllValues);

    % Graph-frequency cutoff estimated from training subjects only
    [idx_low,idx_high] = GSPTB_GetPSDMedianCutoff(X_OI,Data.FD_TimeCourses{id_REST1},...
        H_Pop,idPop,id_session,Cfg.Data.ScrubbingThreshold);

    % Rectangular graph filters
    G.Low = GSPTB_CreateRectangularFilter(idx_low,nHarmonics);
    G.High = GSPTB_CreateRectangularFilter(idx_high,nHarmonics);

    % Population filtering operators
    P_Pop = zeros(nHarmonics,nHarmonics,numel(Filters));

    for id_filter = 1:numel(Filters)
        P_Pop(:,:,id_filter) = H_Pop.U*G.(Filters{id_filter})*H_Pop.U';
    end

    % Initializes feature arrays
    %
    % Each array contains the information for all subjects and, where
    % applicable, both filter types. For functional connectivity we only need
    % the off-diagonal upper triangle of the matrix, while for second-order
    % moment quantities we also retain the diagonal elements
    Cx_All = nan(nEdgesUpper,S,'single');
    FC_All = nan(nEdges,S,'single');

    DC_All = nan(nEdgesUpper,S,numel(Filters),'single');
    EC_All = nan(nEdgesUpper,S,numel(Filters),'single');
    CT_All = nan(nEdgesUpper,S,numel(Filters),'single');

    FC_Pop_All = nan(nEdges,S,numel(Filters),'single');
    FC_Ind_All = nan(nEdges,S,numel(Filters),'single');

    Delta_FC_All = nan(nEdges,S,numel(Filters),'single');
    FC_Error_All = nan(nEdges,S,numel(Filters),'single');

    % Computes features across subjects
    for s = 1:S

        s

        % Subject-specific harmonics
        H_Ind = GSPTB_GetPopulationConnectomeHarmonics(Data.Adj{s});

        % Resting-state activity
        X = double(X_OI{s}{id_session}');

        % Scrubbing
        FD = Data.FD_TimeCourses{id_REST1}{s};
        X = X(:,FD < Cfg.Data.ScrubbingThreshold);

        % Unfiltered features

        % Second-order moment
        Cx = (X*X')/size(X,2);

        % Functional connectivity
        FC = corr(X');

        % Vectorized representations
        FC_All(:,s) = single(jUpperTriMatToVec(FC));
        Cx_All(:,s) = single([diag(Cx);jUpperTriMatToVec(Cx)]);

        % Filter-dependent features
        for id_filter = 1:numel(Filters)

            % Subject-specific filtering operator
            P_Ind = H_Ind.U*G.(Filters{id_filter})*H_Ind.U';

            % Filtered activity
            Y_Pop = P_Pop(:,:,id_filter) * X;
            Y_Ind = P_Ind * X;

            % Functional error
            E = Y_Pop - Y_Ind;

            % Number of retained time points
            T = size(Y_Ind,2);

            % Second-order moment decomposition

            % Difference between population- and individual-filtered
            % second-order moments
            DC = (Y_Pop*Y_Pop')/T - (Y_Ind*Y_Ind')/T;

            % Cross-term
            CT = (Y_Ind*E' + E*Y_Ind')/T;

            % Second-order moment of the functional error
            EC = (E*E')/T;

            % Vectorized representations
            DC_All(:,s,id_filter) = single([diag(DC);jUpperTriMatToVec(DC)]);
            CT_All(:,s,id_filter) = single([diag(CT);jUpperTriMatToVec(CT)]);
            EC_All(:,s,id_filter) = single([diag(EC);jUpperTriMatToVec(EC)]);

            % Functional connectivity features
            FC_Pop = corr(Y_Pop');
            FC_Ind = corr(Y_Ind');
            FC_Error = corr(E');
            Delta_FC = FC_Pop - FC_Ind;

            FC_Pop_All(:,s,id_filter) = single(jUpperTriMatToVec(FC_Pop));
            FC_Ind_All(:,s,id_filter) = single(jUpperTriMatToVec(FC_Ind));
            Delta_FC_All(:,s,id_filter) = single(jUpperTriMatToVec(Delta_FC));
            FC_Error_All(:,s,id_filter) = single(jUpperTriMatToVec(FC_Error));
        end
    end

    % Collects feature matrices
    Features = struct;
    Features.FC = FC_All;
    Features.Cx = Cx_All;
    Features.DC_Low = squeeze(DC_All(:,:,1));
    Features.DC_High = squeeze(DC_All(:,:,2));
    Features.EC_Low = squeeze(EC_All(:,:,1));
    Features.EC_High = squeeze(EC_All(:,:,2));
    Features.CT_Low = squeeze(CT_All(:,:,1));
    Features.CT_High = squeeze(CT_All(:,:,2));
    Features.FC_Pop_Low = squeeze(FC_Pop_All(:,:,1));
    Features.FC_Pop_High = squeeze(FC_Pop_All(:,:,2));
    Features.FC_Ind_Low = squeeze(FC_Ind_All(:,:,1));
    Features.FC_Ind_High = squeeze(FC_Ind_All(:,:,2));
    Features.Delta_FC_Low = squeeze(Delta_FC_All(:,:,1));
    Features.Delta_FC_High = squeeze(Delta_FC_All(:,:,2));
    Features.Error_FC_Low = squeeze(FC_Error_All(:,:,1));
    Features.Error_FC_High = squeeze(FC_Error_All(:,:,2));

    % Prediction parameters
    yTrain = Data.FactorScores.Cognition(idPop);
    yTest = Data.FactorScores.Cognition(idTest);

    % Family structure restricted to outer training subjects
    FamilyKernel_Train = Data.FamilyKernel(idPop,idPop);

    % Test-set nuisance covariate
    CovTest = Data.MeanFD_NonScrubbedOnly(idTest,id_REST1);

    % Family-aware inner folds
    InnerSplits = GSPTB_CreateInnerPredictionSplits(FamilyKernel_Train,PredOpts);

    % Kernel construction and kernel ridge regression across quantities
    for m = 1:numel(FeatureNames)

        % Current feature of interest
        feature = FeatureNames{m};

        % Train and test kernels
        [K_Train,K_Test] = ...
            GSPTB_Create_TrainTestKernels(Features.(feature),idPop,idTest);

        % Kernel ridge regression
        KRR_Results{k}.(feature) = GSPTB_KernelRidgeRegression_TrainTest(...
            K_Train,K_Test,yTrain,yTest,...
            CovTest,InnerSplits,PredOpts);
    end

    % Stores fold identities
    KRR_Results{k}.TrainSubjectID = idPop;
    KRR_Results{k}.TestSubjectID = idTest;

    % Clears fold-specific feature arrays
    clear Features
    clear FC_All Cx_All
    clear DC_All EC_All CT_All
    clear FC_Pop_All FC_Ind_All Delta_FC_All FC_Error_All
    clear A_Pop H_Pop P_Pop G
end



%% 2. Prediction performance
%
% Collects out-of-fold predictions across the 10 outer folds and computes
% prediction performance for each of the 16 feature sets.
%
% Two complementary R2 estimates are retained:
%
%   Fold-wise R2:
%       Prediction performance evaluated separately within each outer
%       test fold.
%
%   Whole-population R2:
%       Prediction performance evaluated after combining the out-of-fold
%       predictions across all subjects.
%
% Each subject contributes exactly one out-of-fold prediction to the
% whole-population estimate
disp('2. Prediction performance...');


% Initializes outputs
nFeatures = numel(FeatureNames);
nFolds = PredSplits.NumFolds;

% Will contain the R2 values obtained from kernel ridge regression, with
% one value per feature and per fold
R2_Distribution = nan(nFeatures,nFolds);

% Will contain the R2 values obtained when predicted values for each
% subject (each subject has a predicted value in only one fold) are
% compared to actual values (whole-population prediction)
R2_WholePop = nan(nFeatures,1);

% Stores complete out-of-fold prediction vectors
PredictionResults = struct;

% Same process across prediction features
for m = 1:nFeatures

    % Current feature
    feature = FeatureNames{m};

    % Will contain one out-of-fold prediction per subject
    yPred_All = nan(S,1);

    % Collects results across outer folds
    for k = 1:nFolds

        % Prediction results for current feature and fold
        tmp = KRR_Results{k}.(feature);

        % Subjects held out in this fold
        idTest = KRR_Results{k}.TestSubjectID;

        % Stores held-out predictions in original subject order
        yPred_All(idTest) = tmp.y_pred;

        % Fold-wise prediction performance
        R2_Distribution(m,k) = tmp.R2;
    end

    % Ensures that every subject received exactly one prediction
    assert(all(~isnan(yPred_All)), ...
        ['Missing out-of-fold predictions for ',feature,'.']);

    % Whole-population prediction performance
    R2_WholePop(m) = GSPTB_ComputeR2(yPred_All,Data.FactorScores.Cognition);

    % Retains complete out-of-fold prediction vector
    PredictionResults.(feature).yPred = yPred_All;
end

% Prediction performance table

% One row per prediction feature
PredictionPerformance = table(string(FeatureNames(:)),R2_WholePop, ...
    'VariableNames',{'Feature','R2_WholePop'});

% Adds fold-wise R2 values
for k = 1:nFolds
    PredictionPerformance.(sprintf('R2_Fold%d',k)) = R2_Distribution(:,k);
end



%% 3. Population versus individual prediction performance
%
% Compares cognition prediction based on FC obtained after filtering with
% population-derived versus subject-specific connectome harmonics.
%
% Analyses are performed separately for low- and high-pass filtering.
%
% Population-versus-individual effects are defined as:
%
%   DeltaR2 = R2_Pop - R2_Ind
%
% such that positive values indicate better prediction using the
% population-derived representation.
%
% Statistical significance is assessed using family-aware permutations.
% Population/individual prediction labels are exchanged jointly for all
% subjects belonging to the same family.
%
% For the frequency-by-representation interaction:
%
%   I = DeltaR2_High - DeltaR2_Low
%
% the same population/individual exchange is applied simultaneously to
% the low- and high-pass predictions
disp('3. Population versus individual prediction performance...');


% Out-of-fold predictions
y = Data.FactorScores.Cognition(:);

yPred_Pop_Low = PredictionResults.FC_Pop_Low.yPred;
yPred_Ind_Low = PredictionResults.FC_Ind_Low.yPred;
yPred_Pop_High = PredictionResults.FC_Pop_High.yPred;
yPred_Ind_High = PredictionResults.FC_Ind_High.yPred;

% Observed prediction performance
R2_Pop_Low = GSPTB_ComputeR2(yPred_Pop_Low,y);
R2_Ind_Low = GSPTB_ComputeR2(yPred_Ind_Low,y);
R2_Pop_High = GSPTB_ComputeR2(yPred_Pop_High,y);
R2_Ind_High = GSPTB_ComputeR2(yPred_Ind_High,y);

% Population minus individual prediction performance
DeltaR2_Low = R2_Pop_Low - R2_Ind_Low;
DeltaR2_High = R2_Pop_High - R2_Ind_High;

% Frequency-by-representation interaction
InteractionR2 = DeltaR2_High - DeltaR2_Low;

% Validation against PredictionPerformance
R2_Pop_Low_Table = PredictionPerformance.R2_WholePop( ...
    PredictionPerformance.Feature == "FC_Pop_Low");

R2_Ind_Low_Table = PredictionPerformance.R2_WholePop( ...
    PredictionPerformance.Feature == "FC_Ind_Low");

R2_Pop_High_Table = PredictionPerformance.R2_WholePop( ...
    PredictionPerformance.Feature == "FC_Pop_High");

R2_Ind_High_Table = PredictionPerformance.R2_WholePop( ...
    PredictionPerformance.Feature == "FC_Ind_High");

assert(abs(R2_Pop_Low-R2_Pop_Low_Table) < 1e-10);
assert(abs(R2_Ind_Low-R2_Ind_Low_Table) < 1e-10);
assert(abs(R2_Pop_High-R2_Pop_High_Table) < 1e-10);
assert(abs(R2_Ind_High-R2_Ind_High_Table) < 1e-10);

% Family structure for permutation testing

% Connected components of the family kernel define exchangeability units
Gfam = graph(Data.FamilyKernel);
FamilyID = conncomp(Gfam);
Families = unique(FamilyID);

% Ensures reproducibility of permutation inference
rng(1);

% Permutation configuration
nPerm = 10000;

DeltaR2_perm_Low = nan(nPerm,1);
DeltaR2_perm_High = nan(nPerm,1);
InteractionR2_perm = nan(nPerm,1);

% Family-aware permutations
for p = 1:nPerm

    % Determines which families undergo a population/individual exchange
    swap = false(S,1);

    for f = 1:numel(Families)
        ids = find(FamilyID == Families(f));

        if rand > 0.5
            swap(ids) = true;
        end
    end

    % Copies of prediction vectors
    ypPL = yPred_Pop_Low;
    ypIL = yPred_Ind_Low;

    ypPH = yPred_Pop_High;
    ypIH = yPred_Ind_High;

    % Population/individual exchange
    %
    % The same exchange is applied simultaneously at low and high
    % frequencies, preserving the paired frequency structure
    tmp = ypPL(swap);
    ypPL(swap) = ypIL(swap);
    ypIL(swap) = tmp;

    tmp = ypPH(swap);
    ypPH(swap) = ypIH(swap);
    ypIH(swap) = tmp;

    % Permuted population-versus-individual effects
    DeltaR2_perm_Low(p) = GSPTB_ComputeR2(ypPL,y) - GSPTB_ComputeR2(ypIL,y);
    DeltaR2_perm_High(p) = GSPTB_ComputeR2(ypPH,y) - GSPTB_ComputeR2(ypIH,y);

    % Permuted interaction
    InteractionR2_perm(p) = DeltaR2_perm_High(p) - DeltaR2_perm_Low(p);
end

% Permutation p-values

% Directional tests:
% H1 = population-derived representation predicts better than
%      subject-specific representation
p_Low = (1 + sum(DeltaR2_perm_Low >= DeltaR2_Low))/(nPerm + 1);

p_High = (1 + sum(DeltaR2_perm_High >= DeltaR2_High))/(nPerm + 1);


% Two-sided frequency-by-representation interaction
p_Interaction = (1 + sum(abs(InteractionR2_perm) >= abs(InteractionR2)))/(nPerm + 1);

% Prediction comparison table
PredictionComparison = table(["Low";"High";"Frequency_x_Representation"],...
    [R2_Pop_Low;R2_Pop_High;NaN], ...
    [R2_Ind_Low;R2_Ind_High;NaN], ...
    [DeltaR2_Low;DeltaR2_High;InteractionR2], ...
    [p_Low;p_High;p_Interaction], ...
    'VariableNames',{'Comparison',...
    'R2_Pop',...
    'R2_Ind',...
    'Effect',...
    'pValue'});



%% 4. Saving all Practical Utility results
% Saves prediction splits, prediction outputs and statistical comparisons
disp('4. Saving all Practical Utility results...');

save(fullfile(ResultsDir,'PracticalUtility_Part3_Results.mat'), ...
    'PredSplits', ...
    'FeatureNames', ...
    'KRR_Results', ...
    'R2_Distribution', ...
    'R2_WholePop', ...
    'PredictionResults', ...
    'PredictionPerformance', ...
    'PredictionComparison', ...
    '-v7.3');