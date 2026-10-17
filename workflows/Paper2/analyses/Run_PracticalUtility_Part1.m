%% 1. Computation of Practical Utility measures
%
% Run time: around 10 minutes
disp('Starting to run Practical Utility analyses (part 1)...');

disp('1. Computation of Practical Utility measures...');

% Number of outer splits
NumOuter_Run = 10;

% Topology and Weight construction schemes to consider
Topology = {'AllEdges','ConsensusEdges'};
Weights = {'AllValues','NZValues'};

% Our analyses consider either low-pass or high-pass graph filtering,
% according to a median power spectral density cutoff
Filters = {'Low','High'};

% Indices required to sample the REST1 LR scan from the data array
id_REST1 = 8;
id_session = 1;

% Loads resting-state data for all subjects; note that we only want
% cortical regions included to match our other analyses
% At this stage, X_OI contains the data for all subjects across both the LR
% and RL acquisitions, respectively indices 1 and 2
X_OI = GSPTB_LoadTaskData(fullfile(DataDir,'Inputs'),'REST1', ...
    Cfg.Atlas.CorticalIndices,Cfg.Atlas.IsCortexOnly);

% Number of harmonics/brain regions
N_Harmonics = Cfg.Atlas.NumRegionsAnalysis;

clear G

% Outer loop
for n_outer = 1:NumOuter_Run

    % Subjects held out from population-harmonic construction
    idTest = PopSplits.Outer(n_outer).Test;

    % Population used to construct the population representation
    idPop = PopSplits.Outer(n_outer).Train;

    % Population connectomes for all construction schemes
    A_Pop = GSPTB_GetPopulationAdjacencyMatrices(Data,Cfg,idPop);

    % Subject-specific harmonics, eigenvalues and Laplacians are obtained
    clear H_Subject

    for s = 1:numel(idTest)
        H_Subject(s) = GSPTB_GetPopulationConnectomeHarmonics(Data.Adj{idTest(s)});
    end

    % Same analysis for all connectome-construction schemes
    for id_topology = 1:numel(Topology)
        for id_weights = 1:numel(Weights)

            % Currently considered scheme details
            topology = Topology{id_topology};
            weights = Weights{id_weights};

            % Population-level harmonics
            H_Pop = GSPTB_GetPopulationConnectomeHarmonics(A_Pop.(topology).(weights));

            % Population-derived graph-frequency cutoff
            % GSPTB_GetPSDMedianCutoff selects only the data from session 
            % id_session (LR acquisition) for all subjects, while the FD
            % time courses in Data are already the LR acquisition ones
            % Computations are done only considering non-scrubbed time
            % points at FD threshold = 0.5 mm
            [idx_low,idx_high] = GSPTB_GetPSDMedianCutoff(X_OI,...
                Data.FD_TimeCourses{id_REST1},H_Pop,idPop, ...
                id_session,Cfg.Data.ScrubbingThreshold);

            % Low- and high-pass graph filters are constructed accordingly
            G.Low = GSPTB_CreateRectangularFilter(idx_low,N_Harmonics);
            G.High = GSPTB_CreateRectangularFilter(idx_high,N_Harmonics);

            % Population filtering operators are built for both filters
            P_Pop_All = zeros(N_Harmonics,N_Harmonics,numel(Filters));

            for id_filter = 1:numel(Filters)
                P_Pop_All(:,:,id_filter) = ...
                    H_Pop.U * G.(Filters{id_filter}) * H_Pop.U';
            end

            % Initializes manuscript-relevant quantities

            % Main quantities of interest
            StructuralMismatch = zeros(numel(idTest),numel(Filters));
            Error_Functional = zeros(numel(idTest),numel(Filters));
            Rel_DeltaC = zeros(numel(idTest),numel(Filters));
            Error_FC = zeros(numel(idTest),numel(Filters));

            % Supplementary quantities (operator and activity-energy
            % quantities)
            StructuralMismatch_Op = zeros(numel(idTest),numel(Filters));
            Energy_Error = zeros(numel(idTest),numel(Filters));
            Energy_Ind = zeros(numel(idTest),numel(Filters));
            Energy_Pop = zeros(numel(idTest),numel(Filters));

            % Supplementary quantities (second-order moment error
            % decomposition)
            Rel_CrossTerm = zeros(numel(idTest),numel(Filters));
            Rel_ErrorSecondMoment = zeros(numel(idTest),numel(Filters));
            Cos_Theta = zeros(numel(idTest),numel(Filters));
            Kappa = zeros(numel(idTest),numel(Filters));

            % Same computations for each held-out subject
            for s = 1:numel(idTest)

                % Gets the harmonics, eigenvalues and Laplacian for the
                % subject
                H_Ind = H_Subject(s);

                % Subject functional data
                X = double(X_OI{idTest(s)}{id_session}');

                % Scrubbing
                FD = Data.FD_TimeCourses{id_REST1}{idTest(s)};
                X = X(:,FD < Cfg.Data.ScrubbingThreshold);

                % Second-order moment matrix of the unfiltered activity
                Cx = (X*X')/size(X,2);

                % Same operations for each filter type
                for id_filter = 1:numel(Filters)

                    % Population- and subject-specific filtering operators
                    P_Pop = P_Pop_All(:,:,id_filter);
                    P_Ind = H_Ind.U * G.(Filters{id_filter}) * H_Ind.U';

                    % Difference between filtering operators
                    D = P_Pop - P_Ind;

                    % Relative structural mismatch with respect to
                    % individual filtering operator
                    % This is termed "Structural mismatch" in the manuscript
                    % and is our first main Practical Utility quantity
                    StructuralMismatch(s,id_filter) = norm(D,'fro')/norm(P_Ind,'fro');

                    StructuralMismatch_Op(s,id_filter) = norm(D,2);

                    % Filtered functional activity using either the
                    % population-level harmonics, or the subject's
                    % individual harmonics
                    Y_Pop = P_Pop * X;
                    Y_Ind = P_Ind * X;

                    % Functional error between both
                    E = Y_Pop - Y_Ind;

                    % Relative functional error with respect to
                    % individually filtered data
                    % This is termed "Functional error" in the manuscript
                    % and is our second main Practical Utility quantity
                    Error_Functional(s,id_filter) = norm(E,'fro')/max(norm(Y_Ind,'fro'),eps);

                    % Verifies that E = D X
                    E_check = D*X;

                    Check_E = norm(E-E_check,'fro')/max(norm(E,'fro'),eps);

                    assert(Check_E < 1e-10, ...
                        'Operator-error identity failed.');

                    % Activity energies, which are supplementary quantities
                    % in the manuscript

                    % Number of time points at hand
                    T = size(X,2);

                    Energy_Pop(s,id_filter) = norm(Y_Pop,'fro')^2/T;
                    Energy_Ind(s,id_filter) = norm(Y_Ind,'fro')^2/T;
                    Energy_Error(s,id_filter) = norm(E,'fro')^2/T;

                    % Verifies the theoretical expression for error energy
                    Error_Energy_Theory = trace(D*Cx*D');

                    Check_Theory = abs(Energy_Error(s,id_filter) - ...
                        Error_Energy_Theory)/max(abs(Energy_Error(s,id_filter)),eps);

                    assert(Check_Theory < 1e-10, ...
                        'Energy-error identity failed.');
                    
                    % Second-order moment error is computed below
                    
                    % Second-order moment matrices for the individually
                    % filtered and population-wise filtered data
                    CI = (Y_Ind*Y_Ind')/T;
                    CP = (Y_Pop*Y_Pop')/T;

                    % Difference between population- and individual-basis
                    % second-order moment matrices
                    DC = CP - CI;

                    % Decomposition of DC into its two constituting terms
                    % CT = cross-term, EC = error second-order moment
                    CT = (Y_Ind*E' + E*Y_Ind')/T;
                    EC = (E*E')/T;

                    % Verifies that DC = CT + EC
                    Check_DeltaC = norm(DC-(CT+EC),'fro')/max(norm(DC,'fro'),eps);

                    assert(Check_DeltaC < 1e-10, ...
                        'Second-order moment decomposition failed.');

                    % To ensure we always divide by something non-zero
                    denCI = max(norm(CI,'fro'),eps);

                    % Relative difference between second-order moment
                    % matrices with respect to the individual-basis
                    % second-order moment matrix
                    % This is termed "Second-order moment error" in the 
                    % manuscript and is our third main Practical Utility 
                    % quantity
                    Rel_DeltaC(s,id_filter) = norm(DC,'fro')/denCI;

                    % The next quantities are further supplementary ones in
                    % terms of the second-order moment error decomposition
                    Rel_CrossTerm(s,id_filter) = norm(CT,'fro')/denCI;

                    Rel_ErrorSecondMoment(s,id_filter) = norm(EC,'fro')/denCI;

                    % Alignment between decomposition terms
                    Cos_Theta(s,id_filter) = trace(CT'*EC)/...
                        max(norm(CT,'fro')*norm(EC,'fro'),eps);

                    % Cancellation index
                    Kappa(s,id_filter) = -2*trace(CT'*EC)/...
                        max(norm(CT,'fro')^2 + norm(EC,'fro')^2,eps);

                    
                    % Functional connectivity error is computed below
                    
                    % Functional connectivity matrices obtained from
                    % individual-filtered and population-filtered
                    % functional data
                    FC_Pop = corr(Y_Pop');
                    FC_Ind = corr(Y_Ind');

                    % Vectorized versions
                    vFC_Pop = jUpperTriMatToVec(FC_Pop);
                    vFC_Ind = jUpperTriMatToVec(FC_Ind);

                    % Relative difference between individual-filtered and
                    % population-filtered functional connectivity vectors
                    % with respect to the individual-filtered data
                    % This is termed "FC error" in the manuscript and is 
                    % our fourth main Practical Utility quantity
                    Error_FC(s,id_filter) = ...
                        norm(vFC_Pop-vFC_Ind) / ...
                        max(norm(vFC_Ind),eps);
                end
            end

            % Stores all manuscript-relevant quantities into P
            P = struct;

            % Main quantities
            P.StructuralMismatch = StructuralMismatch;
            P.Error_Functional = Error_Functional;
            P.Rel_DeltaC = Rel_DeltaC;
            P.Error_FC = Error_FC;

            % Supplementary quantities
            P.StructuralMismatch_Op = StructuralMismatch_Op;

            P.Energy_Error = Energy_Error;
            P.Energy_Ind = Energy_Ind;
            P.Energy_Pop = Energy_Pop;

            P.Rel_CrossTerm = Rel_CrossTerm;
            P.Rel_ErrorSecondMoment = Rel_ErrorSecondMoment;
            
            P.Cos_Theta = Cos_Theta;
            P.Kappa = Kappa;

            % Frequency-cutoff information
            P.Cutoff.Low = idx_low;
            P.Cutoff.High = idx_high;
            P.CutoffIndex = idx_high(1);

            % We performed the computations for a given outer split and
            % connectome construction scheme
            PracticalUtility.Outer(n_outer).(topology).(weights) = P;
        end
    end

    % Subject identities associated with this outer split
    PracticalUtility.Outer(n_outer).SubjectID = idTest;
end

% Filters
PracticalUtility.Filters = Filters;



%% 2. Construction of Practical Utility data tables
% (a) Observation-level table across outer splits, subjects, construction
% schemes and filtering types
% (b) Subject-level table where repeated observations have been averaged
% across outer splits
%
% Table (b) is only used for visualizations, but statistical analyses are
% always performed from Table (a) to preserve the Outer x Subject
% observations
disp('2. Construction of Practical Utility data tables...');

% Quantities retained for manuscript analyses
MetricNames = {'StructuralMismatch','Error_Functional','Rel_DeltaC','Error_FC', ...
    'StructuralMismatch_Op','Energy_Error','Energy_Ind','Energy_Pop', ...
    'Rel_CrossTerm','Rel_ErrorSecondMoment','Cos_Theta','Kappa'};

% The lines below build the observation-level table (a)

% Canonical table
% One row per Outer x Subject x Topology x Weights x Filter
PracticalUtilityStats = table;

% Going across outer splits...
for o = 1:NumOuter_Run

    % Subjects left-out in that split
    Subjects = PracticalUtility.Outer(o).SubjectID;

    % Goes across connectome construction schemes
    for id_topology = 1:numel(Topology)

        topology = Topology{id_topology};

        for id_weights = 1:numel(Weights)

            weights = Weights{id_weights};

            % Gets the data we need
            P = PracticalUtility.Outer(o).(topology).(weights);

            % Same process across filter types
            for id_filter = 1:numel(Filters)

                % Which filtering type
                filter = Filters{id_filter};

                % How many subjects
                nSubjects = numel(Subjects);

                % Identifiers
                Outer_col = repmat(o,nSubjects,1);
                Subject_col = Subjects(:);
                Topology_col = repmat(string(topology),nSubjects,1);
                Weights_col = repmat(string(weights),nSubjects,1);
                Filter_col = repmat(string(filter),nSubjects,1);

                % Manuscript-relevant quantities

                % Main quantities
                StructuralMismatch_col = P.StructuralMismatch(:,id_filter);
                Error_Functional_col = P.Error_Functional(:,id_filter);
                Rel_DeltaC_col = P.Rel_DeltaC(:,id_filter);
                Error_FC_col = P.Error_FC(:,id_filter);

                % Supplementary quantities
                StructuralMismatch_Op_col = P.StructuralMismatch_Op(:,id_filter);

                Energy_Error_col = P.Energy_Error(:,id_filter);
                Energy_Ind_col = P.Energy_Ind(:,id_filter);
                Energy_Pop_col = P.Energy_Pop(:,id_filter);

                Rel_CrossTerm_col = P.Rel_CrossTerm(:,id_filter);
                Rel_ErrorSecondMoment_col = P.Rel_ErrorSecondMoment(:,id_filter);

                Cos_Theta_col = P.Cos_Theta(:,id_filter);
                Kappa_col = P.Kappa(:,id_filter);

                % Adds current observations
                T = table(Outer_col,Subject_col,Topology_col,Weights_col,Filter_col,...
                    StructuralMismatch_col,Error_Functional_col,Rel_DeltaC_col,Error_FC_col,...
                    StructuralMismatch_Op_col, ...
                    Energy_Error_col,Energy_Ind_col,Energy_Pop_col, ...
                    Rel_CrossTerm_col,Rel_ErrorSecondMoment_col, ...
                    Cos_Theta_col,Kappa_col, ...
                    'VariableNames',{'Outer','Subject','Topology','Weights','Filter',...
                    'StructuralMismatch','Error_Functional','Rel_DeltaC','Error_FC', ...
                    'StructuralMismatch_Op', ...
                    'Energy_Error','Energy_Ind','Energy_Pop', ...
                    'Rel_CrossTerm','Rel_ErrorSecondMoment', ...
                    'Cos_Theta','Kappa'});

                PracticalUtilityStats = ...
                    [PracticalUtilityStats; T];
            end
        end
    end
end

% Checks that all expected observations have been retained
ExpectedRows = 0;

for o = 1:NumOuter_Run
    ExpectedRows = ExpectedRows + ...
        numel(PracticalUtility.Outer(o).SubjectID)*numel(Topology) * ...
        numel(Weights) * numel(Filters);
end

assert(height(PracticalUtilityStats) == ExpectedRows, ...
    'Unexpected number of rows in PracticalUtilityStats.');

% Categorical grouping variables
PracticalUtilityStats.Outer = categorical(PracticalUtilityStats.Outer);
PracticalUtilityStats.Topology = categorical(PracticalUtilityStats.Topology);
PracticalUtilityStats.Weights = categorical(PracticalUtilityStats.Weights);
PracticalUtilityStats.Filter = categorical(PracticalUtilityStats.Filter);

% The lines below build the specific table (b)

% Converts manuscript-relevant quantities to long format
PULong = stack(PracticalUtilityStats,MetricNames, ...
    'NewDataVariableName','Value','IndexVariableName','Metric');

% Each subject can occur in more than one held-out outer split. Values are
% therefore averaged across occurrences to obtain one value per subject,
% construction scheme, filter and metric
GroupVars_Subject = {'Subject','Topology','Weights','Filter','Metric'};
[G,TGroups] = findgroups(PULong(:,GroupVars_Subject));
NumObservations = splitapply(@numel,PULong.Value,G);
MeanAcrossOuter = splitapply(@mean,PULong.Value,G);

% Builds the resulting table
PracticalUtilitySubjectMeans = [TGroups,table(NumObservations,MeanAcrossOuter, ...
    'VariableNames',{'NumObservations','MeanAcrossOuter'})];



%% 3. Statistical analysis of frequency-dependent differences
disp('3. Statistical analysis of frequency-dependent differences...');

% High - Low differences are assessed for all manuscript-relevant
% Practical Utility quantities

% Separates low- and high-pass observations
TLow = PULong(PULong.Filter == "Low",:);
THigh = PULong(PULong.Filter == "High",:);

% Guarantees exact pairing between filtering conditions
SortVars = {'Outer','Subject','Topology','Weights','Metric'};

TLow = sortrows(TLow,SortVars);
THigh = sortrows(THigh,SortVars);

assert(height(TLow) == height(THigh));
assert(isequal(TLow.Outer,THigh.Outer));
assert(isequal(TLow.Subject,THigh.Subject));
assert(isequal(TLow.Topology,THigh.Topology));
assert(isequal(TLow.Weights,THigh.Weights));
assert(isequal(TLow.Metric,THigh.Metric));

% Positive Delta indicates larger values for high-pass filtering
Delta = THigh.Value - TLow.Value;

% Stores paired frequency-dependent differences
PracticalUtilityFilterDifferences = table(TLow.Outer,TLow.Subject,...
    TLow.Topology,TLow.Weights,...
    TLow.Metric,TLow.Value,THigh.Value,Delta,...
    'VariableNames',{'Outer','Subject','Topology','Weights',...
    'Metric','Low','High','Delta'});

% Mixed-model analysis of frequency-dependent differences is conducted
% below

% Will contain mixed-model results for all construction schemes and metrics
PracticalUtilityFilterContrasts = table;

% Same process across construction schemes
for id_topology = 1:numel(Topology)

    topology = Topology{id_topology};

    for id_weights = 1:numel(Weights)

        weights = Weights{id_weights};

        % Same process across all metrics
        for id_metric = 1:numel(MetricNames)

            metric = MetricNames{id_metric};

            % Selects the current construction scheme and metric
            % Retains entries across outer splits and subjects in terms of
            % Delta (high - low)
            T = PracticalUtilityFilterDifferences( ...
                PracticalUtilityFilterDifferences.Topology == topology & ...
                PracticalUtilityFilterDifferences.Weights == weights & ...
                PracticalUtilityFilterDifferences.Metric == metric,:);

            % Statistical table
            TDelta = table(categorical(T.Outer),categorical(T.Subject), ...
                T.Delta, ...
                'VariableNames',{'Outer','Subject','Delta'});

            % Fits:
            % Delta_os = beta_0 + u_o + w_s + epsilon_os
            mdl = fitlme(TDelta, ...
                'Delta ~ 1 + (1|Outer) + (1|Subject)');

            % Fixed-effect outputs for beta_0
            C = mdl.Coefficients;

            Estimate = C.Estimate(1);
            SE = C.SE(1);
            DF = C.DF(1);
            tStat = C.tStat(1);
            pValue = C.pValue(1);
            CI_Lower = C.Lower(1);
            CI_Upper = C.Upper(1);

            % Empirical summaries of paired High - Low differences
            Median_Delta = median(T.Delta);
            P025_Delta = prctile(T.Delta,2.5);
            P975_Delta = prctile(T.Delta,97.5);

            % Variance components
            [Psi,MSE] = covarianceParameters(mdl);

            Var_Outer = Psi{1};
            Var_Subject = Psi{2};
            Var_Residual = MSE;

            % Standardized model estimate
            SD_Total = sqrt(Var_Outer + Var_Subject + Var_Residual);
            d_LME = Estimate/SD_Total;

            % Tangibility is assessed only for the four main normalized
            % Practical Utility quantities
            if ismember(metric,{'StructuralMismatch','Error_Functional', ...
                    'Rel_DeltaC','Error_FC'})
            
                TangibilityThreshold = Cfg.PracticalUtility.TangibilityThreshold;
            
                if abs(Estimate) >= TangibilityThreshold
                    Tangible = "Yes";
                else
                    Tangible = "No";
                end
            
            else
                TangibilityThreshold = NaN;
                Tangible = "Not assessed";
            end

            % Final summary row
            NewRow = table(string(topology),string(weights),string(metric), ...
                Median_Delta,P025_Delta,P975_Delta, ...
                Estimate,CI_Lower,CI_Upper,SE,DF,tStat,pValue, ...
                Var_Outer,Var_Subject,Var_Residual, ...
                SD_Total,d_LME,TangibilityThreshold,Tangible, ...
                'VariableNames',{'Topology','Weights','Metric', ...
                'Median_Delta','P025_Delta','P975_Delta', ...
                'Estimate','CI_Lower','CI_Upper','SE','DF','tStat','pValue', ...
                'Var_Outer','Var_Subject','Var_Residual', ...
                'SD_Total','d_LME','TangibilityThreshold','Tangible'});

            PracticalUtilityFilterContrasts = ...
                [PracticalUtilityFilterContrasts; NewRow];
        end
    end
end



%% 4. Saving all Practical Utility results
% Saves all outputs required for downstream analyses and figure generation
disp('4. Saving all Practical Utility results...');

save(fullfile(ResultsDir,'PracticalUtility_Part1_Results.mat'), ...
    'PracticalUtility', ...
    'PracticalUtilityStats', ...
    'PracticalUtilitySubjectMeans', ...
    'PracticalUtilityFilterDifferences', ...
    'PracticalUtilityFilterContrasts', ...
    '-v7.3');


