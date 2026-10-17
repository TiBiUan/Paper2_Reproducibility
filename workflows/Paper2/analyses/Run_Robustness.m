%% 1. Computation of measures of Robustness
% (a) Computation of all possible Robustness measures in the devised
% data splitting scheme (across outer splits, Monte Carlo realizations,
% Topology/Weight construction choices and sample sizes)
% (b) Generation of a table containing information about median and "worst
% percentile" (P5 or P95) for all measures, across outer splits, Monte
% Carlo realizations, Topology/Weight construction choices and sample sizes
%
% Running time on full data: around 20 minutes
disp('Starting to run Robustness analyses...');

disp('1. Computation of measures of Robustness...');

% Will contain our summary data
RobustnessStats = table;

% Number of connectome harmonics (R=216 in the paper)
N_harmonics = Cfg.Atlas.NumRegionsAnalysis;

% Investigated sub-cases regarding the construction of population-level
% connectomes
Topology = {'AllEdges','ConsensusEdges'};
Weights = {'AllValues','NZValues'};

% Below is the core computational scheme

% Loops across outer splits...
for n_outer = 1:PopSplits.NumOuterSplits

    % For each outer split, we focus only on training subjects, and sample
    % subsets n_MC times (Monte Carlo realizations)
    for n_MC = 1:PopSplits.NumMonteCarlo

        % For each Monte Carlo realization, we build population connectome
        % representations using increasingly large subsets of subjects in
        % nested fashion
        for ss = 1:numel(PopSplits.SampleSizes)

            % These are the indices of the subjects to consider to create
            % the two population bases (Set A and Set B)
            idA = PopSplits.Outer(n_outer).MonteCarlo(n_MC).SampleA{ss};
            idB = PopSplits.Outer(n_outer).MonteCarlo(n_MC).SampleB{ss};

            % idB is constructed so that it is empty if one cannot find
            % enough remaining training samples to build it (as the samples
            % must be fully disjoint across sets); in such cases, we just
            % move to the next candidate sample size without further
            % calculations
            if isempty(idB)
                continue;
            end

            % For both sets, computes all candidate connectomes according
            % to the examined topology and weight approaches
            A_Pop_SampleA = GSPTB_GetPopulationAdjacencyMatrices(Data,Cfg,idA);
            A_Pop_SampleB = GSPTB_GetPopulationAdjacencyMatrices(Data,Cfg,idB);

            % Performs similar computations across all these choices
            for id_topology = 1:numel(Topology)
                for id_weights = 1:numel(Weights)

                    % Samples the currently considered construction choices
                    topology = Topology{id_topology};
                    weights = Weights{id_weights};

                    % H1 and H2 contain the Laplacian operator, harmonics
                    % and associated graph frequencies for both sets of
                    % subjects
                    H1 = GSPTB_GetPopulationConnectomeHarmonics(...
                        A_Pop_SampleA.(topology).(weights));
                    
                    H2 = GSPTB_GetPopulationConnectomeHarmonics(...
                        A_Pop_SampleB.(topology).(weights));

                    % Contains all evaluations for the current outer
                    % split/Monte Carlo realization/sample
                    % size/construction scheme
                    Results = struct;

                    % The lines below perform the Robustness computations
                    % per se

                    % Harmonic-wise correspondence
                    Results.HarmonicSimilarity.Index = ...
                        GSPTB_ComputeHarmonicSimilarity(H1.U,H2.U,'index');
                    
                    % Because matching is directional, we match Set B 
                    % harmonics to Set A and Set A harmonics to Set B, 
                    % then average the two results; the same symmetrization 
                    % is used for all directional metrics below
                    Results.HarmonicSimilarity.Match = ...
                        (GSPTB_ComputeHarmonicSimilarity(H1.U,H2.U,'match')+...
                        GSPTB_ComputeHarmonicSimilarity(H2.U,H1.U,'match'))/2;

                    % Spectral transfer matrices, used to compute several
                    % quantities of interest
                    T12 = GSPTB_ComputeSpectralTransfer(H1.U,H2.U);
                    T21 = T12';

                    % Spectral bandwidth using indices (directional, so averaged)
                    BW_AB = GSPTB_ComputeSpectralTransfer_SpectralBandwidth(...
                        T12,H1.Lambda,H2.Lambda,'index');
                    
                    BW_BA = GSPTB_ComputeSpectralTransfer_SpectralBandwidth(...
                        T21,H2.Lambda,H1.Lambda,'index');

                    Results.SpectralTransfer.BandwidthIndex = (BW_AB+BW_BA)/2;
                    
                    % Leakage entropy (directional, so averaged)
                    LE_AB = GSPTB_ComputeSpectralTransfer_LeakageEntropy(T12);
                    LE_BA = GSPTB_ComputeSpectralTransfer_LeakageEntropy(T21);

                    Results.SpectralTransfer.LeakageEntropy = (LE_AB + LE_BA)/2;

                    % Projection error (directional, so subsequent required
                    % window is averaged)
                    PE_AB = GSPTB_ComputeProjectionError(T12,...
                        Cfg.Population.ProjectionWindows);
                    
                    PE_BA = GSPTB_ComputeProjectionError(T21,...
                        Cfg.Population.ProjectionWindows);

                    % Required windows (window sizes enabling to reach a 
                    % projection error lower than a threshold [0.1 in the paper])
                    % All half-window sizes from 0 to 215 are considered
                    RW_AB = GSPTB_ComputeProjectionError_RequiredSpectralWindow(...
                        PE_AB,Cfg.Population.ProjectionWindows,...
                        Cfg.Population.ProjectionTolerance);
                    
                    RW_BA = GSPTB_ComputeProjectionError_RequiredSpectralWindow(...
                        PE_BA,Cfg.Population.ProjectionWindows,...
                        Cfg.Population.ProjectionTolerance);
                    
                    % We may keep NaN values if no window size allows the
                    % criterion to be satisfied, so we isolate valid 
                    % projection windows
                    ValidProjectionWindow = ~isnan(RW_AB) & ~isnan(RW_BA);

                    Results.ValidProjectionWindow = ValidProjectionWindow;
                    Results.RequiredProjectionWindow = (RW_AB + RW_BA)/2;

                    % Subspace similarity (non-directional)
                    % Half-window sizes from 0 to 150 (36 values) are assessed
                    [~,Results.RequiredSubspaceWindow] = ...
                        GSPTB_ComputeSubspaceSimilarityAcrossWindows(...
                        H1.U,H2.U,...
                        Cfg.Population.SubspaceWindows,...
                        Cfg.Population.SubspaceSimilarityThreshold);

                    % Again, we keep track of valid windows as well
                    Results.ValidSubspaceWindow = ~isnan(Results.RequiredSubspaceWindow);

                    % Laplacian reconstruction error (directional, so
                    % averaged)
                    Rec_AB = GSPTB_ComputeLaplacianReconstructionError(H1.U,H2.L);
                    Rec_BA = GSPTB_ComputeLaplacianReconstructionError(H2.U,H1.L);
                    
                    Results.LaplacianReconstructionError = (Rec_AB + Rec_BA)/2;

                    % Saving sample size for the considered occurrence
                    Results.N = PopSplits.SampleSizes(ss);

                    % All the computed measures are saved in the
                    % "Robustness" structure
                    Robustness.Outer(n_outer).MonteCarlo(n_MC)...
                        .SampleSize(ss).(topology).(weights) = Results;
                    
                    % The lines below fill in the summary table

                    % For window-based results, we need to make sure that
                    % we only consider non-NaN entries to compute medians
                    % and percentiles
                    RPW = Results.RequiredProjectionWindow;
                    RPW = RPW(~isnan(RPW));
                    
                    RSW = Results.RequiredSubspaceWindow;
                    RSW = RSW(~isnan(RSW));
                    
                    if isempty(RPW)
                        Median_ProjErr_ReqWin = NaN;
                        P95_ProjErr_ReqWin = NaN;
                    else
                        Median_ProjErr_ReqWin = median(RPW);
                        P95_ProjErr_ReqWin = prctile(RPW,95);
                    end
                    
                    if isempty(RSW)
                        Median_SubSim_ReqWin = NaN;
                        P95_SubSim_ReqWin = NaN;
                    else
                        Median_SubSim_ReqWin = median(RSW);
                        P95_SubSim_ReqWin = prctile(RSW,95);
                    end

                    % Adds the resulting values to our statistics table
                    NewRow = table(n_outer,n_MC,Results.N,...
                        string(topology),string(weights), ...
                        ...
                        median(Results.HarmonicSimilarity.Index),...
                        prctile(Results.HarmonicSimilarity.Index,5),...
                        ...
                        median(Results.HarmonicSimilarity.Match),...
                        prctile(Results.HarmonicSimilarity.Match,5),...
                        sum(Results.HarmonicSimilarity.Match>=0.9)/N_harmonics,...
                        sum(Results.HarmonicSimilarity.Match>=0.95)/N_harmonics,...
                        ...
                        median(Results.SpectralTransfer.BandwidthIndex),...
                        prctile(Results.SpectralTransfer.BandwidthIndex,95),...
                        ...
                        median(Results.SpectralTransfer.LeakageEntropy),...
                        prctile(Results.SpectralTransfer.LeakageEntropy,95),...
                        ...
                        Median_ProjErr_ReqWin,...
                        P95_ProjErr_ReqWin,...
                        sum(Results.ValidProjectionWindow)/N_harmonics,...
                        ...
                        Median_SubSim_ReqWin,...
                        P95_SubSim_ReqWin,...
                        sum(Results.ValidSubspaceWindow)/N_harmonics,...
                        ...
                        Results.LaplacianReconstructionError,...
                        ...
                        'VariableNames',{'Outer','MC','N',...
                        'Topology','Weights', ...
                        'Median_HarmSim_Index','P05_HarmSim_Index',...
                        'Median_HarmSim_Match','P05_HarmSim_Match',...
                        'Frac_HarmSim_Match_090','Frac_HarmSim_Match_095',...
                        'Median_SpecTran_BanInd','P95_SpecTran_BanInd',...
                        'Median_SpecTran_LeakEnt','P95_SpecTran_LeakEnt',...
                        'Median_ProjErr_ReqWin','P95_ProjErr_ReqWin',...
                        'FracValid_ProjErr_ReqWin',...
                        'Median_SubSim_ReqWin','P95_SubSim_ReqWin',...
                        'FracValid_SubSim_ReqWin',...
                        'LaplacianReconstructionError'});

                    RobustnessStats = [RobustnessStats; NewRow];
                 end
            end
        end
    end

    % Checkpoint save, saving the output after each outer split
    save('Robustness_Checkpoint.mat', ...
    'Robustness','RobustnessStats','n_outer','PopSplits','-v7.3');
end

% Makes the required values categorical
RobustnessStats.Outer = categorical(RobustnessStats.Outer);
RobustnessStats.MC = categorical(RobustnessStats.MC);
RobustnessStats.Topology = categorical(RobustnessStats.Topology);
RobustnessStats.Weights = categorical(RobustnessStats.Weights);

% Creates a variable reflective of the interaction between outer splits and
% Monte Carlo realizations
RobustnessStats.Trajectory = categorical( ...
    "O" + string(RobustnessStats.Outer) + ...
    "_MC" + string(RobustnessStats.MC));



%% 2. Computation of descriptive statistics
% Summarizes each Robustness measure across outer splits and Monte Carlo
% realizations, separately for each sample size and connectome construction
% scheme
%
% For each case, computes the mean, standard deviation, median, 2.5th and
% 97.5th percentiles, and number of valid observations
%
% Required inputs for this subsection: RobustnessStats
disp('2. Computation of descriptive statistics...');

% Metrics to summarize
MetricNames = { ...
    'Median_HarmSim_Index', ...
    'P05_HarmSim_Index', ...
    'Median_HarmSim_Match', ...
    'P05_HarmSim_Match', ...
    'Frac_HarmSim_Match_090', ...
    'Frac_HarmSim_Match_095', ...
    'Median_SpecTran_BanInd', ...
    'P95_SpecTran_BanInd', ...
    'Median_SpecTran_LeakEnt', ...
    'P95_SpecTran_LeakEnt', ...
    'Median_ProjErr_ReqWin', ...
    'P95_ProjErr_ReqWin', ...
    'FracValid_ProjErr_ReqWin', ...
    'Median_SubSim_ReqWin', ...
    'P95_SubSim_ReqWin', ...
    'FracValid_SubSim_ReqWin', ...
    'LaplacianReconstructionError'};

NValues = unique(RobustnessStats.N);
TopologyValues = unique(RobustnessStats.Topology);
WeightValues = unique(RobustnessStats.Weights);

% Final long-format table
RobustnessDescriptives = table;

for id_topology = 1:numel(TopologyValues)

    topology = TopologyValues(id_topology);

    for id_weights = 1:numel(WeightValues)

        weights = WeightValues(id_weights);

        for id_N = 1:numel(NValues)

            N = NValues(id_N);

            % Current subset:
            % all Outer x MC realizations for this N and construction scheme
            idx = RobustnessStats.Topology == topology & ...
                RobustnessStats.Weights == weights & ...
                RobustnessStats.N == N;

            for id_metric = 1:numel(MetricNames)

                metric = MetricNames{id_metric};

                x = RobustnessStats.(metric)(idx);
                x = x(~isnan(x));

                if isempty(x)
                    Mean = NaN;
                    SD = NaN;
                    Median = NaN;
                    P025 = NaN;
                    P975 = NaN;
                    NumObservations = 0;
                else
                    Mean = mean(x);
                    SD = std(x);
                    Median = median(x);
                    P025 = prctile(x,2.5);
                    P975 = prctile(x,97.5);
                    NumObservations = numel(x);
                end

                NewRow = table( ...
                    string(metric), ...
                    N, ...
                    topology, ...
                    weights, ...
                    NumObservations, ...
                    Mean, ...
                    SD, ...
                    Median, ...
                    P025, ...
                    P975, ...
                    'VariableNames',{ ...
                    'Metric', ...
                    'N', ...
                    'Topology', ...
                    'Weights', ...
                    'NumObservations', ...
                    'Mean', ...
                    'SD', ...
                    'Median', ...
                    'P025', ...
                    'P975'});

                RobustnessDescriptives = ...
                    [RobustnessDescriptives; NewRow];
            end
        end
    end
end



%% 3. Robustness contrast analyses at literature-inferred values
% Here, we consider specific contrasts to examine further, justified by the
% literature
%
% Required inputs for this subsection: RobustnessStats, Cfg
%
% Only computes the outputs for one structural connectome construction
% scheme
disp('3. Robustness contrast analyses at literature-inferred values...');

% Reference population-harmonic construction scheme, for which we sample
% the obtained values
T = RobustnessStats(RobustnessStats.Topology == "AllEdges" & ...
    RobustnessStats.Weights == "AllValues",:);

% Sample-size contrasts motivated by the literature
N_Contrasts = Cfg.Population.RobustnessContrasts;

% Main Robustness measures retained for the sample size analyses
% reported in the manuscript
Metrics = {'Median_HarmSim_Match', ...
    'P95_ProjErr_ReqWin', ...
    'LaplacianReconstructionError'};

% Direction of improvement:
% +1 = larger values indicate better robustness
% -1 = smaller values indicate better robustness
ImprovementDirection = [1;   % Median_HarmSim_Match
    -1;  % P95_ProjErr_ReqWin
    -1]; % LaplacianReconstructionError

% Will contain final model outputs
RobustnessContrasts = table;

% Optionally retain the actual delta values and fitted models
RobustnessContrastData = struct;
RobustnessContrastModels = struct;

% Initializes the extra array that will contain Laplacian error results
RelativeLaplacianReconstructionError = struct;

% We apply the same approach for each metric of interest...
for id_metric = 1:numel(Metrics)

    metric = Metrics{id_metric};
    direction = ImprovementDirection(id_metric);
    TangibilityThreshold = Cfg.Robustness.Tangibility.(metric);

    % ... and for each contrast to consider
    for id_contrast = 1:size(N_Contrasts,1)

        % N values to look at
        N1 = N_Contrasts(id_contrast,1);
        N2 = N_Contrasts(id_contrast,2);

        % Extract both sample size-associated value sets
        T1 = T(T.N == N1,:);
        T2 = T(T.N == N2,:);

        % Ensure identical ordering across paired trajectories
        T1 = sortrows(T1,{'Outer','MC'});
        T2 = sortrows(T2,{'Outer','MC'});

        % Basic integrity checks
        assert(height(T1) == height(T2));
        assert(isequal(T1.Outer,T2.Outer));
        assert(isequal(T1.MC,T2.MC));

        % Positive Delta always means improvement
        if direction == 1
            Delta = T2.(metric) - T1.(metric);
        else
            Delta = T1.(metric) - T2.(metric);
        end

        % For Laplacian reconstruction error, we also compute the outcomes
        % for relative values with respect to the N1
        if strcmp(metric,'LaplacianReconstructionError')

            % Relative improvement with respect to the "starting point" [%]
            DeltaRel = 100 * (T1.LaplacianReconstructionError-...
                T2.LaplacianReconstructionError)./T1.LaplacianReconstructionError;
            
            % Delta-level statistical table
            TDeltaRel = table(categorical(T1.Outer),categorical(T1.MC),...
            DeltaRel,'VariableNames',{'Outer','MC','DeltaRel'});

            % Fits the model of interest to the data
            mdl_rel = fitlme(TDeltaRel,'DeltaRel ~ 1 + (1|Outer)');

            % Gets beta coefficient and its standard error, as well as
            % confidence interval and p-value (fixed effect outcome for beta_0)
            tmp = mdl_rel.Coefficients;
            RelativeLaplacianReconstructionError(id_contrast).Est = tmp.Estimate;
            RelativeLaplacianReconstructionError(id_contrast).SE = tmp.SE;
            RelativeLaplacianReconstructionError(id_contrast).CI = ...
                [tmp.Lower,tmp.Upper];
            RelativeLaplacianReconstructionError(id_contrast).pValue = tmp.pValue;
            
            % Gets random effect information (variance) and residual
            % variance
            [PsiRel,MSERel] = covarianceParameters(mdl_rel);
            RelativeLaplacianReconstructionError(id_contrast).VarOuter = PsiRel{1};
            RelativeLaplacianReconstructionError(id_contrast).VarResidual = MSERel;

            % Keeps track of the model per se
            RelativeLaplacianReconstructionError(id_contrast).Model = mdl_rel;
        end

        % Delta-level statistical table
        TDelta = table(categorical(T1.Outer),categorical(T1.MC),Delta, ...
            'VariableNames',{'Outer','MC','Delta'});

        % Fits Delta_om = beta_0 + w_o + epsilon_om
        mdl = fitlme(TDelta,'Delta ~ 1 + (1|Outer)');

        % Fixed effect outputs (i.e., for beta_0)
        C = mdl.Coefficients;

        % Selects specific information (estimate, standard error, degrees
        % of freedom, t-statistic, p-value, confidence interval
        Estimate = C.Estimate(1);
        SE = C.SE(1);
        DF = C.DF(1);
        tStat = C.tStat(1);
        pValue = C.pValue(1);
        CI_Lower = C.Lower(1);
        CI_Upper = C.Upper(1);

        % Variance components
        [Psi,MSE] = covarianceParameters(mdl);

        % One random effect only (Outer) has its variance retained, and we
        % also keep track of residual variance
        Var_Outer = Psi{1};
        Var_Residual = MSE;

        SD_Outer = sqrt(Var_Outer);
        SD_Residual = sqrt(Var_Residual);

        % Intraclass correlation
        ICC_Outer = Var_Outer/(Var_Outer+Var_Residual);

        % Empirical Delta summaries
        Mean_Delta = mean(Delta);
        SD_Delta = std(Delta);
        Median_Delta = median(Delta);
        P025_Delta = prctile(Delta,2.5);
        P975_Delta = prctile(Delta,97.5);

        % Determines whether the model-estimated improvement reaches the
        % predefined metric-specific threshold for a tangible change
        if strcmp(metric,'LaplacianReconstructionError')
            Tangible = RelativeLaplacianReconstructionError(id_contrast).Est >= ...
                TangibilityThreshold;
        else
            Tangible = Estimate >= TangibilityThreshold;
        end

        % How often did the individual trajectory improve?
        FracPositive = mean(Delta>0);
        FracNonNegative = mean(Delta>=0);

        % Identifier for current contrast
        Contrast = string(N1) + "_to_" + string(N2);
        ContrastField = "N" + string(N1) + "_to_" + string(N2);

        % Final summary row
        NewRow = table( ...
            string(metric), ...
            N1, ...
            N2, ...
            Contrast, ...
            Mean_Delta, ...
            SD_Delta, ...
            Median_Delta, ...
            P025_Delta, ...
            P975_Delta, ...
            FracPositive, ...
            FracNonNegative, ...
            Estimate, ...
            SE, ...
            DF, ...
            tStat, ...
            pValue, ...
            CI_Lower, ...
            CI_Upper, ...
            Var_Outer, ...
            Var_Residual, ...
            SD_Outer, ...
            SD_Residual, ...
            ICC_Outer, ...
            Tangible,...
            'VariableNames',{ ...
            'Metric', ...
            'N1', ...
            'N2', ...
            'Contrast', ...
            'Mean_Delta', ...
            'SD_Delta', ...
            'Median_Delta', ...
            'P025_Delta', ...
            'P975_Delta', ...
            'FracPositive', ...
            'FracNonNegative', ...
            'Estimate', ...
            'SE', ...
            'DF', ...
            'tStat', ...
            'pValue', ...
            'CI_Lower', ...
            'CI_Upper', ...
            'Var_Outer', ...
            'Var_Residual', ...
            'SD_Outer', ...
            'SD_Residual', ...
            'ICC_Outer',...
            'Tangible'});

        RobustnessContrasts = ...
            [RobustnessContrasts; NewRow];

        % Retains complete Delta data
        RobustnessContrastData.(metric).(ContrastField).Delta = Delta;
        RobustnessContrastData.(metric).(ContrastField).Outer = TDelta.Outer;
        RobustnessContrastData.(metric).(ContrastField).MC = TDelta.MC;

        % Retains fitted model
        RobustnessContrastModels.(metric).(ContrastField) = mdl;
    end
end

% Added for effect size computations 
RobustnessContrasts.StdTotal = sqrt(RobustnessContrasts.Var_Outer + ...
    RobustnessContrasts.Var_Residual);

RobustnessContrasts.d_LME = RobustnessContrasts.Estimate ./ ...
    RobustnessContrasts.StdTotal;



%% 4. Analyses of convergence rates
% Quantifies the marginal gain in Robustness per additional subject and
% fits power-law convergence models to each {Outer,MC} trajectory for the
% three Robustness measures retained in the manuscript
%
% Required inputs for this subsection: Metrics, ImprovementDirection,
% NValues, RobustnessStats, Cfg
%
% Only computes the outputs for the reference structural connectome
% construction scheme (AllEdges-AllValues)
disp('4. Analyses of convergence rates...');

% Data for the construction scheme of interest
T = RobustnessStats(RobustnessStats.Topology == "AllEdges" & ...
    RobustnessStats.Weights == "AllValues",:);

% Reference ordering of {Outer,MC} trajectories
TRef = T(T.N == NValues(1),:);
TRef = sortrows(TRef,{'Outer','MC'});

NumTrajectories = height(TRef);

% Checks that identical trajectories are available at every sample size
for n = 1:numel(NValues)

    Tn = T(T.N == NValues(n),:);
    Tn = sortrows(Tn,{'Outer','MC'});

    assert(height(Tn) == NumTrajectories, ...
        'Sample sizes do not contain identical numbers of trajectories.');

    assert(isequal(Tn.Outer,TRef.Outer) && isequal(Tn.MC,TRef.MC), ...
        'Sample sizes do not contain identical Outer x MC trajectories.');
end

% Middle value between each sample size pair
Nmid = (NValues(1:end-1) + NValues(2:end))/2;

% Difference in sample size between consecutive values
dN = diff(NValues);

% Y will contain metric values for all {o,mc} realizations across sample
% sizes
% n_metrics x trajectory (Outer x MC) x number of sample sizes
Y = nan(numel(Metrics),NumTrajectories,numel(NValues));

% The lines below perform the computation of gain rates (Supplementary Figure 1A)

% Same analysis for each metric
for met = 1:length(Metrics)  
    
    % Gets the values for the metric of interest across sample sizes
    for n = 1:numel(NValues)

        Tn = T(T.N == NValues(n),:);
        Tn = sortrows(Tn,{'Outer','MC'});
    
        Y(met,:,n) = Tn.(Metrics{met});
    end

    % Computes gain rate for all realizations as well as summarizing
    % statistics across them
    % Gain rate is expressed as change in the metric per additional subject
    GainRate{met} = ImprovementDirection(met) * diff(squeeze(Y(met,:,:)),1,2)./dN';
    GainMed{met}  = median(GainRate{met},1);
    GainP025{met} = prctile(GainRate{met},2.5,1);
    GainP975{met} = prctile(GainRate{met},97.5,1);
end

% The lines below perform convergence fitting (Figure 2B, Supplementary Figure 1B)

% Smallest studied population size
N0 = Cfg.Population.RobustnessSampleSizes(1);

% Same process for each metric
for met = 1:length(Metrics)

    % The model is adjusted depending on which direction an improvement
    % occurs for, but has the same parameters in both cases
    if ImprovementDirection(met)==1
        ft = fittype( ...
            'Yinf - a*x^(-p)', ...
            'independent','x', ...
            'coefficients',{'Yinf','a','p'});
    else
        ft = fittype( ...
            'Yinf + a*x^(-p)', ...
            'independent','x', ...
            'coefficients',{'Yinf','a','p'});
    end    

    % Same process for each {o,mc} trajectory across sample sizes
    for traj = 1:size(Y,2)
        
        % We get x and y data (respectively sample size values and
        % associated metric values)
        x = NValues(:);
        y = squeeze(Y(met,traj,:));

        opts = fitoptions(ft);
        opts.Lower = [0 0 0];

        % If we consider median harmonic similarity, then we cannot exceed
        % 1, so we set the upper limit accordingly
        if strcmp(Metrics{met},'Median_HarmSim_Match')
            opts.Upper = [1 Inf Inf];
            opts.StartPoint = [min(1, max(y)),max(y)-min(y),0.5];

        % Else, we have no clear assumption regarding upper boundaries
        else
            opts.Upper = [Inf Inf Inf];
            opts.StartPoint = [max(0, min(y)),max(y)-min(y),0.5];
        end

        % Fitting process per se
        [cf,gof] = fit(x,y,ft,opts);

        % Extracting parameters
        PowerLawParams.Yinf(met,traj) = cf.Yinf;
        PowerLawParams.a(met,traj) = cf.a;
        PowerLawParams.p(met,traj) = cf.p;

        % Extracting fitting quality measures 
        FittingQuality.R2(met,traj) = gof.rsquare;
        FittingQuality.AR2(met,traj) = gof.adjrsquare;
        FittingQuality.RMSE(met,traj) = gof.rmse;

        % Population size at which the fitted distance to the asymptote is
        % 10% of its value at the smallest investigated sample size (N0)
        PowerLawParams.N90(met,traj) = N0*10.^(1./cf.p);
    end
end



%% 5. Saving all Robustness results
% Saves all outputs required for downstream analyses and figure generation
disp('5. Saving all Robustness results...');

save(fullfile(ResultsDir,'Robustness_Results.mat'), ...
    'Robustness', ...
    'RobustnessStats', ...
    'RobustnessDescriptives', ...
    'RobustnessContrasts', ...
    'RobustnessContrastData', ...
    'RobustnessContrastModels', ...
    'RelativeLaplacianReconstructionError', ...
    'GainRate', ...
    'GainMed', ...
    'GainP025', ...
    'GainP975', ...
    'Nmid', ...
    'PowerLawParams', ...
    'FittingQuality', ...
    '-v7.3');