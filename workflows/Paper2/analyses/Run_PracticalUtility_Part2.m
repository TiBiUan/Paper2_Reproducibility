%% 1. Computation of Practical Utility across population sample sizes
%
% Practical Utility is evaluated within the same population splitting scheme
% used for the other pillars:
%
% Outer split
%     -> Monte Carlo realization
%         -> Population sample size
%             -> Held-out subjects
%
% The following quantities are retained:
%
% Main manuscript:
%   StructuralMismatch
%   Error_Functional
%   Rel_DeltaC
%   Error_FC
%
% Supplementary analyses:
%   Rel_CrossTerm
%   Rel_ErrorSecondMoment
%   Cos_Theta
%   Kappa
%
% The current analyses use only the reference population construction:
%   Topology = AllEdges
%   Weights  = AllValues
%
% Run time: approximately 50 minutes
disp('Starting to run Practical Utility analyses (part 2)...');

disp('1. Computation of Practical Utility across population sample sizes...');

% The lines below set the configuration

% Number of outer splits and Monte Carlo realizations
NumOuter_Run = PopSplits.NumOuterSplits;
NumMC_Run = PopSplits.NumMonteCarlo;

% In this part, we only consider our main connectome construction scheme
Topology = {'AllEdges'};
Weights = {'AllValues'};

% We consider low- and high-frequency graph filtering
Filters = {'Low','High'};

% Index of REST1 session in the loaded HCP functional data structure
id_REST1 = 8;

% First resting-state run
id_session = 1;

% Each cell corresponds to one subject
% The second cell index corresponds to the recording/run
X_OI = GSPTB_LoadTaskData(fullfile(DataDir,'Inputs'),'REST1',...
    Cfg.Atlas.CorticalIndices,Cfg.Atlas.IsCortexOnly);

% Number of investigated sample sizes
nSS = numel(PopSplits.SampleSizes);

% Number of harmonics
nHarmonics = size(Data.Adj{1},1);

% The lines below initialize the output structure
PracticalUtility_SS.SampleSizes = PopSplits.SampleSizes;
PracticalUtility_SS.Filters = Filters;
PracticalUtility_SS.Topology = Topology;
PracticalUtility_SS.Weights = Weights;

% Makes sure we do not conserve the "G" grouping variable from before,
% because we use G for filters below
clear G

% Across outer splits...
for n_outer = 1:NumOuter_Run

    % Subjects held out in this outer split
    idTest = PopSplits.Outer(n_outer).Test;

    % Stores subject identities once
    PracticalUtility_SS.Outer(n_outer).SubjectID = idTest;

    % Subject-specific quantities

    % Structural harmonics, eigenvalues and Laplacian for held-out subjects
    clear H_Subject

    for s = 1:numel(idTest)
        H_Subject(s) = GSPTB_GetPopulationConnectomeHarmonics(Data.Adj{idTest(s)});
    end

    % Functional covariance quantities that do not depend on N are
    % initialized
    X_Subject = cell(numel(idTest),1);
    Cx_Subject = cell(numel(idTest),1);
    Ccorr_Subject = cell(numel(idTest),1);

    % Same process for each held-out subject
    for s = 1:numel(idTest)

        % Functional time series for the subject at hand
        X = double(X_OI{idTest(s)}{id_session}');

        % Scrubbing of the data
        FD = Data.FD_TimeCourses{id_REST1}{idTest(s)};
        X = X(:,FD < Cfg.Data.ScrubbingThreshold);
        X_Subject{s} = X;

        % Covariance matrix for the subject at hand
        Cx_Subject{s} = (X*X')/size(X,2);

        % Correlation-derived covariance for FC estimation
        Xc = X - mean(X,2);
        Ccorr_Subject{s} = Xc*Xc';
    end

    % Across Monte Carlo realizations...
    for n_MC = 1:NumMC_Run

        n_MC

        % Across sample sizes...
        for ss = 1:nSS

            % Subjects used to construct population connectome
            idPop = PopSplits.Outer(n_outer).MonteCarlo(n_MC).SampleA{ss};

            % Population adjacency matrices
            A_Pop = GSPTB_GetPopulationAdjacencyMatrices(Data,Cfg,idPop);

            % Population harmonics, eigenvalues and Laplacian for our
            % construction scheme of interest
            H_Pop = GSPTB_GetPopulationConnectomeHarmonics(A_Pop.AllEdges.AllValues);

            % Frequency cutoffs and filters
            [idx_low,idx_high] = GSPTB_GetPSDMedianCutoff(X_OI,...
                Data.FD_TimeCourses{id_REST1},H_Pop,idPop,id_session,...
                Cfg.Data.ScrubbingThreshold);

            G.Low = GSPTB_CreateRectangularFilter(idx_low,nHarmonics);
            G.High = GSPTB_CreateRectangularFilter(idx_high,nHarmonics);

            % Population filtering operators
            P_Pop_All = zeros(nHarmonics,nHarmonics,numel(Filters));

            for id_filter = 1:numel(Filters)
                P_Pop_All(:,:,id_filter) = H_Pop.U*G.(Filters{id_filter})*H_Pop.U';
            end

            % Initializes the outputs

            % Number of left-out subjects
            nSubjects = numel(idTest);

            % All four main quantities
            StructuralMismatch = zeros(nSubjects,numel(Filters));
            Error_Functional = zeros(nSubjects,numel(Filters));
            Rel_DeltaC = zeros(nSubjects,numel(Filters));
            Error_FC = zeros(nSubjects,numel(Filters));

            % Supplementary quantities
            Rel_CrossTerm = zeros(nSubjects,numel(Filters));
            Rel_ErrorSecondMoment = zeros(nSubjects,numel(Filters));
            Cos_Theta = zeros(nSubjects,numel(Filters));
            Kappa = zeros(nSubjects,numel(Filters));

            % Held-out subjects
            for s = 1:nSubjects

                % Harmonics, eigenvalues and Laplacian for the subject at
                % hand
                H_Ind = H_Subject(s);

                % Second-order moment and FC
                Cx = Cx_Subject{s};
                Ccorr = Ccorr_Subject{s};

                % Same for all filter types...
                for id_filter = 1:numel(Filters)

                    % Population and individual filtering operators
                    P_Pop = squeeze(P_Pop_All(:,:,id_filter));
                    P_Ind = H_Ind.U*G.(Filters{id_filter})*H_Ind.U';

                    % Difference between operators
                    D = P_Pop - P_Ind;

                    % Structural mismatch
                    StructuralMismatch(s,id_filter) = norm(D,'fro')/norm(P_Ind,'fro');

                    % Functional propagation
                    CI = P_Ind * Cx * P_Ind';
                    CP = P_Pop * Cx * P_Pop';
                    EC = D * Cx * D';
                    CT = P_Ind * Cx * D' + D * Cx * P_Ind';

                    % Difference in second-order moment
                    DC = CP - CI;

                    % To avoid dividing by zero
                    denCI = max(norm(CI,'fro'),eps);

                    % Functional error
                    Error_Functional(s,id_filter) = sqrt(trace(EC)/max(trace(CI),eps));

                    % Relative second-order moment difference
                    Rel_DeltaC(s,id_filter) = norm(DC,'fro')/denCI;

                    % Other quantities
                    Rel_CrossTerm(s,id_filter) = norm(CT,'fro')/denCI;
                    Rel_ErrorSecondMoment(s,id_filter) = norm(EC,'fro')/denCI;
                    Cos_Theta(s,id_filter) = trace(CT'*EC)/...
                        max(norm(CT,'fro')*norm(EC,'fro'),eps);

                    Kappa(s,id_filter) = -2*trace(CT'*EC)/...
                        max(norm(CT,'fro')^2 + norm(EC,'fro')^2,eps);

                    % Functional connectivity error
                    CFC_Pop = P_Pop*Ccorr*P_Pop';
                    CFC_Ind = P_Ind*Ccorr*P_Ind';

                    dPop = sqrt(max(diag(CFC_Pop),eps));
                    FC_Pop = CFC_Pop./(dPop*dPop');

                    dInd = sqrt(max(diag(CFC_Ind),eps));
                    FC_Ind = CFC_Ind./(dInd*dInd');

                    % Vectorized versions
                    vFC_Pop = jUpperTriMatToVec(FC_Pop);
                    vFC_Ind = jUpperTriMatToVec(FC_Ind);

                    % Relative FC error
                    Error_FC(s,id_filter) = norm(vFC_Pop-vFC_Ind)/max(norm(vFC_Ind),eps);
                end
            end

            % Stores the results
            PU = struct;

            PU.StructuralMismatch = StructuralMismatch;
            PU.Error_Functional = Error_Functional;
            PU.Rel_DeltaC = Rel_DeltaC;
            PU.Error_FC = Error_FC;
           
            PU.Rel_CrossTerm = Rel_CrossTerm;
            PU.Rel_ErrorSecondMoment = Rel_ErrorSecondMoment;
            PU.Cos_Theta = Cos_Theta;
            PU.Kappa = Kappa;

            % Saves in PracticalUtility_SS
            PracticalUtility_SS.Outer(n_outer).MonteCarlo(n_MC).SampleSize(ss)...
                .AllEdges.AllValues = PU;
        end
    end

    % Saves outer split independently
    PU_Outer = PracticalUtility_SS.Outer(n_outer);

    save(sprintf('PracticalUtility_SS_Outer%02d.mat',n_outer),'PU_Outer','-v7.3');
end



%% 2. Construction of Practical Utility data tables
% (a) Observation-level table across outer splits, Monte Carlo
% realizations, population sample sizes, subjects and filtering types
% (b) Global descriptive table across outer splits, Monte Carlo
% realizations and subjects for each population sample size and filter
disp('2. Construction of Practical Utility data tables...');

% Quantities retained for manuscript analyses
MetricNames = {'StructuralMismatch','Error_Functional','Rel_DeltaC','Error_FC', ...
    'Rel_CrossTerm','Rel_ErrorSecondMoment','Cos_Theta','Kappa'};

% The lines below build the observation-level table (a)

% Canonical table
% One row per Outer x MC x N x Subject x Filter
PracticalUtilitySSStats = table;

% Goes across outer splits
for o = 1:NumOuter_Run

    % Subjects held out in this outer split
    Subjects = PracticalUtility_SS.Outer(o).SubjectID;
    nSubjects = numel(Subjects);

    % Goes across Monte Carlo realizations
    for mc = 1:NumMC_Run

        % Goes across population sample sizes
        for ss = 1:nSS

            N = PracticalUtility_SS.SampleSizes(ss);

            % Gets the data we need
            P = PracticalUtility_SS.Outer(o).MonteCarlo(mc) ...
                .SampleSize(ss).AllEdges.AllValues;

            % Same process across filtering types
            for id_filter = 1:numel(Filters)

                filter = Filters{id_filter};

                % Identifiers
                Outer_col = repmat(o,nSubjects,1);
                MC_col = repmat(mc,nSubjects,1);
                N_col = repmat(N,nSubjects,1);
                Subject_col = Subjects(:);
                Filter_col = repmat(string(filter),nSubjects,1);

                % Main quantities
                StructuralMismatch_col = ...
                    P.StructuralMismatch(:,id_filter);

                Error_Functional_col = ...
                    P.Error_Functional(:,id_filter);

                Rel_DeltaC_col = ...
                    P.Rel_DeltaC(:,id_filter);

                Error_FC_col = ...
                    P.Error_FC(:,id_filter);

                % Supplementary quantities
                Rel_CrossTerm_col = ...
                    P.Rel_CrossTerm(:,id_filter);

                Rel_ErrorSecondMoment_col = ...
                    P.Rel_ErrorSecondMoment(:,id_filter);

                Cos_Theta_col = ...
                    P.Cos_Theta(:,id_filter);

                Kappa_col = ...
                    P.Kappa(:,id_filter);

                % Adds current observations
                T = table(Outer_col,MC_col,N_col,Subject_col,Filter_col, ...
                    StructuralMismatch_col,Error_Functional_col, ...
                    Rel_DeltaC_col,Error_FC_col, ...
                    Rel_CrossTerm_col,Rel_ErrorSecondMoment_col, ...
                    Cos_Theta_col,Kappa_col, ...
                    'VariableNames',{'Outer','MC','N','Subject','Filter', ...
                    'StructuralMismatch','Error_Functional', ...
                    'Rel_DeltaC','Error_FC', ...
                    'Rel_CrossTerm','Rel_ErrorSecondMoment', ...
                    'Cos_Theta','Kappa'});

                PracticalUtilitySSStats = ...
                    [PracticalUtilitySSStats; T];
            end
        end
    end
end

% Checks that all expected observations have been retained
ExpectedRows = 0;

for o = 1:NumOuter_Run
    ExpectedRows = ExpectedRows + ...
        numel(PracticalUtility_SS.Outer(o).SubjectID) * ...
        NumMC_Run * nSS * numel(Filters);
end

assert(height(PracticalUtilitySSStats) == ExpectedRows, ...
    'Unexpected number of rows in PracticalUtilitySSStats.');

% Categorical grouping variables
PracticalUtilitySSStats.Outer = ...
    categorical(PracticalUtilitySSStats.Outer);

PracticalUtilitySSStats.MC = ...
    categorical(PracticalUtilitySSStats.MC);

PracticalUtilitySSStats.Filter = ...
    categorical(PracticalUtilitySSStats.Filter);

% The lines below build the global descriptive table (b)

% Converts manuscript-relevant quantities to long format
PULong = stack(PracticalUtilitySSStats,MetricNames, ...
    'NewDataVariableName','Value','IndexVariableName','Metric');

% Summarizes values across outer splits, Monte Carlo realizations and
% subjects for each population sample size, filter and metric
GroupVars = {'N','Filter','Metric'};

[G,TGroups] = findgroups(PULong(:,GroupVars));

NumObservations = splitapply(@numel,PULong.Value,G);
Mean = splitapply(@mean,PULong.Value,G);
SD = splitapply(@std,PULong.Value,G);
Median = splitapply(@median,PULong.Value,G);
P025 = splitapply(@(x) prctile(x,2.5),PULong.Value,G);
P975 = splitapply(@(x) prctile(x,97.5),PULong.Value,G);

PracticalUtilitySSDescriptives = [TGroups,table(NumObservations,Mean,SD,Median,P025,P975)];



%% 3. Sample size contrasts
%
% Tests changes in the four primary Practical Utility quantities across
% literature-motivated increases in population sample size.
%
% All four quantities are error measures. Delta is therefore defined as:
%
%   Delta = Value(N1) - Value(N2)
%
% such that positive values indicate improvement (reduced error) with
% increasing population sample size.
%
% Mixed-effects model:
%
%   Delta ~ 1 + (1|Outer) + (1|OuterMC) + (1|Subject)
disp('3. Sample size contrasts...');


%% Configuration

% Sample size contrasts
N_Contrasts = Cfg.Population.PracticalUtilityContrasts;

% Primary Practical Utility quantities
MetricNames_Main = {'StructuralMismatch','Error_Functional','Rel_DeltaC','Error_FC'};

% Tangibility threshold for normalized relative-error measures
TangibilityThreshold = Cfg.PracticalUtility.TangibilityThreshold;

% Initialize results table
PracticalUtilitySSContrasts = table;

% Contrasts
for id_filter = 1:numel(Filters)

    filter = Filters{id_filter};

    % Data linked to the current filtering type
    TFilter = PracticalUtilitySSStats(PracticalUtilitySSStats.Filter == filter,:);

    % Same across metrics...
    for id_metric = 1:numel(MetricNames_Main)

        metric = MetricNames_Main{id_metric};

        % Same for each contrast...
        for id_contrast = 1:size(N_Contrasts,1)

            % Sample sizes of interest
            N1 = N_Contrasts(id_contrast,1);
            N2 = N_Contrasts(id_contrast,2);

            % Extract observations at the two sample sizes
            T1 = TFilter(TFilter.N == N1, ...
                {'Outer','MC','Subject',metric});

            T2 = TFilter(TFilter.N == N2, ...
                {'Outer','MC','Subject',metric});

            % Sorts observations to ensure exact pairing
            PairVars = {'Outer','MC','Subject'};

            T1 = sortrows(T1,PairVars);
            T2 = sortrows(T2,PairVars);

            % Checks exact correspondence between observations
            assert(height(T1) == height(T2), ...
                'Different numbers of observations at N1 and N2.');

            assert(isequal(T1.Outer,T2.Outer) && ...
                isequal(T1.MC,T2.MC) && ...
                isequal(T1.Subject,T2.Subject), ...
                'Observations are not correctly paired across sample sizes.');

            % Computes the change in quantity
            %
            % All quantities are errors:
            % positive Delta = reduction in error = improvement

            Delta = T1.(metric) - T2.(metric);

            % Mixed-effects table
            TDelta = table(T1.Outer,T1.MC,T1.Subject,Delta,...
                'VariableNames',{'Outer','MC','Subject','Delta'});

            % Monte Carlo realizations are nested within outer splits
            TDelta.OuterMC = categorical(string(TDelta.Outer) + "_" + string(TDelta.MC));

            % Mixed-effects model
            mdl = fitlme(TDelta, ...
                'Delta ~ 1 + (1|Outer) + (1|OuterMC) + (1|Subject)');

            % Fixed-effect statistics
            C = mdl.Coefficients;

            Estimate = C.Estimate(1);
            SE = C.SE(1);
            DF = C.DF(1);
            tStat = C.tStat(1);
            pValue = C.pValue(1);
            CI_Lower = C.Lower(1);
            CI_Upper = C.Upper(1);

            % Empirical distribution of paired changes
            MedianDelta = median(Delta);
            P025Delta = prctile(Delta,2.5);
            P975Delta = prctile(Delta,97.5);

            % Tangibility
            Tangible = abs(Estimate) >= TangibilityThreshold;

            if Tangible
                TangibleString = "Yes";
            else
                TangibleString = "No";
            end

            % Stores the results
            NewRow = table(string(filter),string(metric),N1,N2,...
                Estimate,SE,DF,tStat,pValue,CI_Lower,CI_Upper,...
                MedianDelta,P025Delta,P975Delta,...
                TangibilityThreshold,TangibleString,...
                'VariableNames',{'Filter','Metric','N1','N2',...
                'Estimate','SE','DF','tStat','pValue','CI_Lower','CI_Upper',...
                'MedianDelta','P025Delta','P975Delta',...
                'TangibilityThreshold','Tangible'});

            PracticalUtilitySSContrasts = ...
                [PracticalUtilitySSContrasts;NewRow];
        end
    end
end



%% 4. Saving all Practical Utility results
% Saves all outputs required for downstream analyses and figure generation
disp('4. Saving all Practical Utility results...');

save(fullfile(ResultsDir,'PracticalUtility_Part2_Results.mat'), ...
    'PracticalUtility_SS', ...
    'PracticalUtilitySSStats', ...
    'PracticalUtilitySSDescriptives', ...
    'PracticalUtilitySSContrasts', ...
    '-v7.3');