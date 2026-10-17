%% 1. Computation of measures of Representativeness
% (a) Computation of all possible Representativeness measures in the devised
% data splitting scheme (across outer splits, Monte Carlo realizations,
% Topology/Weight construction choices and sample sizes)
%
% Run time for full dataset:
% Harmonic similarity: around 15 hours
% Bandwidth index: around 20 minutes
% Leakage entropy: around 30 minutes
% Subspace similarity: around 20 minutes
disp(['Starting to run Representativeness analyses for ',RepCfg.MetricBaseName,'...']);

disp('1. Computation of measures of Representativeness...');

% Local copies of important parameters for the data splitting scheme
NumOuter_Run = PopSplits.NumOuterSplits;
NumMC_Run = PopSplits.NumMonteCarlo;

% Metric-specific configuration
switch RepCfg.MetricBaseName

    case 'HarmonicSimilarityMatch'

        % Larger values indicate better correspondence
        RepCfg.Central.Name = 'Median';
        RepCfg.Central.Function = @(x) median(x,4);

        % Lower tail reflects the least representative harmonics
        RepCfg.Tail.Percentile = 5;

        % Better Representativeness means larger similarity
        RepCfg.ImprovementDirection = [1 1];

        RepCfg.TangibilityStatistic = 'Central';
        RepCfg.TangibilityThreshold = 0.01;

    % For all the cases below, larger metric values mean worse
    % Representativeness
    case {'BandwidthIndex','LeakageEntropy','ProjectionErrorRequiredWindows', ...
            'LaplacianReconstructionError','SubspaceSimilarityRequiredWindows'}

        % Smaller values indicate better correspondence
        RepCfg.Central.Name = 'Median';
        RepCfg.Central.Function = @(x) median(x,4);

        % Upper tail reflects the least representative harmonics
        RepCfg.Tail.Percentile = 95;

        RepCfg.ImprovementDirection = [-1 -1];

        if strcmp(RepCfg.MetricBaseName,'BandwidthIndex')
            RepCfg.TangibilityStatistic = 'Central';
            RepCfg.TangibilityThreshold = 0.5;
        elseif strcmp(RepCfg.MetricBaseName,'LeakageEntropy')
            RepCfg.TangibilityStatistic = 'Central';
            RepCfg.TangibilityThreshold = 0.01*log(Cfg.Atlas.NumRegionsAnalysis);
        elseif strcmp(RepCfg.MetricBaseName,'SubspaceSimilarityRequiredWindows')
            RepCfg.TangibilityStatistic = 'Tail';
            RepCfg.TangibilityThreshold = 1;
        else
            RepCfg.TangibilityStatistic = '';
            RepCfg.TangibilityThreshold = NaN;
        end

    otherwise

        error('Unknown Representativeness metric: %s', ...
            RepCfg.MetricBaseName);
end

% The lines below add to RepCfg required fields for subsequent table
% generation as well as to keep track of computed statistics

% Tail statistic
RepCfg.Tail.Name = sprintf('P%d',RepCfg.Tail.Percentile);
RepCfg.Tail.Function = @(x) prctile(x,RepCfg.Tail.Percentile,4);

% Automatically generated variable names
RepCfg.Central.VariableName = sprintf('%s_%s',RepCfg.Central.Name,...
    RepCfg.MetricBaseName);

RepCfg.Tail.VariableName = sprintf('%s_%s',RepCfg.Tail.Name,...
    RepCfg.MetricBaseName);

RepCfg.MetricNames = {RepCfg.Central.VariableName,RepCfg.Tail.VariableName};

% Reference population-connectome construction scheme
RepCfg.ReferenceTopology = "AllEdges";
RepCfg.ReferenceWeights = "AllValues";

% Sample size contrasts
RepCfg.N_Contrasts = Cfg.Population.RepresentativenessContrasts;

% Population-connectome construction schemes considered in the paper;
% AllEdges-AllValues is the reference scheme used for the main analyses
Topology = {'AllEdges','ConsensusEdges'};
Weights = {'AllValues','NZValues'};

% Pre-initializes the cell array that will contain our Representativeness
% measure
RepMeas = cell(PopSplits.NumOuterSplits,1);

% Creating properly-sized arrays that will be filled with actual values
% below (single data type)
for o = 1:PopSplits.NumOuterSplits
    nSub = PopSplits.Outer(o).NumTest;

    for it = 1:numel(Topology)
        for iw = 1:numel(Weights)

            topology = Topology{it};
            weights = Weights{iw};
            RepMeas{o}.(topology).(weights) = ...
                nan(PopSplits.NumMonteCarlo,numel(PopSplits.SampleSizes),...
                nSub,Cfg.Atlas.NumRegionsAnalysis,'single');
        end
    end
end

% Keeps track of which are the left-out subjects for each outer split
RepMeasSubjects = cell(PopSplits.NumOuterSplits,1);

for o = 1:PopSplits.NumOuterSplits
    RepMeasSubjects{o} = PopSplits.Outer(o).Test;
end

% First, we have the loop across outer splits
for n_outer = 1:NumOuter_Run

    clear H_Subject

    % We need to keep track of our test subjects this time
    idTest = PopSplits.Outer(n_outer).Test;

    % Handy variables
    nSS = numel(PopSplits.SampleSizes);
    nTest = numel(idTest);
    nTopology = numel(Topology);
    nWeights = numel(Weights);

    % Computes connectome harmonics for all held-out subjects
    for s = 1:numel(idTest)
        H_Subject(s) = GSPTB_GetPopulationConnectomeHarmonics(Data.Adj{idTest(s)});
    end

    % For each outer split, we focus only on training subjects and consider
    % Monte Carlo realizations one after the other
    for n_MC = 1:NumMC_Run

        % For each realization, we build population connectome
        % representations from increasingly large subsets
        for ss = 1:nSS

            % Indices of the subjects used to construct the current
            % population representation
            idPop = PopSplits.Outer(n_outer).MonteCarlo(n_MC).SampleA{ss};
            
            % Gets all four population-connectome representations
            A_Pop = GSPTB_GetPopulationAdjacencyMatrices(Data,Cfg,idPop);

            % Same process for all population-connectome construction schemes
            for id_topology = 1:nTopology
                for id_weights = 1:nWeights

                    topology = Topology{id_topology};
                    weights = Weights{id_weights};

                    % Population-level harmonics and eigenvalues
                    H_Pop = GSPTB_GetPopulationConnectomeHarmonics(...
                        A_Pop.(topology).(weights));

                    % Compare this population representation
                    % with every held-out individual
                    for s = 1:nTest

                        % Harmonics for the presently considered subject
                        H_Ind = H_Subject(s);

                        % Computations are performed for the metric
                        % currently being considered
                        switch RepCfg.MetricBaseName

                            % Matched harmonic-wise similarity
                            case 'HarmonicSimilarityMatch'
                        
                                RepValue = GSPTB_ComputeHarmonicSimilarity( ...
                                    H_Pop.U,H_Ind.U,'match');
                        
                            % Bandwidth in index units
                            case 'BandwidthIndex'
                        
                                STM = GSPTB_ComputeSpectralTransfer(H_Pop.U,H_Ind.U);
                        
                                RepValue = ...
                                    GSPTB_ComputeSpectralTransfer_SpectralBandwidth( ...
                                    STM,H_Pop.Lambda,H_Ind.Lambda,'index');
                        
                            % Leakage entropy
                            case 'LeakageEntropy'
                        
                                STM = GSPTB_ComputeSpectralTransfer(H_Pop.U,H_Ind.U);
                        
                                RepValue = ...
                                    GSPTB_ComputeSpectralTransfer_LeakageEntropy(STM);
                        
                            % Number of windows required for projection
                            % error to become lower than tolerance level
                            case 'ProjectionErrorRequiredWindows'
                                
                                STM = GSPTB_ComputeSpectralTransfer(H_Pop.U,H_Ind.U);
                                PE = GSPTB_ComputeProjectionError( ...
                                    STM,Cfg.Population.ProjectionWindows);

                                RepValue = ...
                                    GSPTB_ComputeProjectionError_RequiredSpectralWindow( ...
                                    PE, ...
                                    Cfg.Population.ProjectionWindows, ...
                                    Cfg.Population.ProjectionTolerance);
                        
                            % Laplacian reconstruction error
                            case 'LaplacianReconstructionError'
                                
                               RepValue = ...
                                    GSPTB_ComputeLaplacianReconstructionError( ...
                                    H_Pop.U,H_Ind.L);
                        
                            % Number of windows required for subspace
                            % similarity to exceed similarity threshold
                            case 'SubspaceSimilarityRequiredWindows'
                                
                                [~,RepValue] = ...
                                    GSPTB_ComputeSubspaceSimilarityAcrossWindows( ...
                                    H_Pop.U,H_Ind.U, ...
                                    Cfg.Population.SubspaceWindows, ...
                                    Cfg.Population.SubspaceSimilarityThreshold);
                            otherwise
                                error('Unknown Representativeness metric: %s', ...
                                    RepCfg.MetricBaseName);
                        end

                        assert(numel(RepValue) == Cfg.Atlas.NumRegionsAnalysis, ...
                            'Representativeness metric must return one value per harmonic.');

                        % RepMeas contains our results for the metric
                        % being currently considered
                        RepMeas{n_outer}.(topology).(weights)(n_MC,ss,s,:) = single(RepValue);
                    end
                end
            end
        end
    end
end



%% 2. Computation of Representativeness descriptive statistics
% (a) Subject-level summaries across harmonics
% For each subject, population realization, sample size and connectome
% construction scheme, computes a central and tail summary of
% Representativeness across harmonics
disp('2. Computation of Representativeness descriptive statistics...');

% Computes how many rows we expect to have in the table
nRows = 0;

for o = 1:NumOuter_Run

    % How many left-out subjects we have in the examined outer split
    nSub = numel(RepMeasSubjects{o});

    % How many rows in total are expected for that outer split
    nRows = nRows + NumMC_Run * numel(PopSplits.SampleSizes) * nSub * ...
        numel(Topology) * numel(Weights);
end

% Preallocation
Outer_col = zeros(nRows,1);
MC_col = zeros(nRows,1);
N_col = zeros(nRows,1);
Subject_col = zeros(nRows,1);
Topology_col = strings(nRows,1);
Weights_col = strings(nRows,1);
Central_col = nan(nRows,1);
Tail_col = nan(nRows,1);

row = 0;

% Builds the table's data
for o = 1:NumOuter_Run

    % Held-out subjects for the current outer split
    Subjects = RepMeasSubjects{o};

    % Same process for all connectome construction schemes
    for it = 1:numel(Topology)
        for iw = 1:numel(Weights)

            % Current connectome construction scheme
            topology = Topology{it};
            weights  = Weights{iw};

            % X has size n_MC x n_samplesizes x n_subjects x n_harmonics
            X = RepMeas{o}.(topology).(weights);
            X = X(1:NumMC_Run,:,:,:);

            % Computes across harmonics (fourth dimension of X)
            Central = RepCfg.Central.Function(X);
            Tail = RepCfg.Tail.Function(X);
            
            % Central and Tail have size n_MC x n_samplesizes x n_subjects
            % after removing the singleton harmonic dimension
            Central = reshape(Central,NumMC_Run,numel(PopSplits.SampleSizes), ...
                numel(Subjects));
            
            Tail = reshape(Tail,NumMC_Run,numel(PopSplits.SampleSizes), ...
                numel(Subjects));

            % We now unwrap the computed values into separate rows
            for mc = 1:NumMC_Run
                for n = 1:numel(PopSplits.SampleSizes)

                    nr = numel(Subjects);

                    % For each outer split, Monte Carlo realization and 
                    % sample size, we always add n_subjects values (how
                    % many are left-out at the considered outer split)
                    idx = row + (1:nr);

                    % Information about splits, realizations, sample sizes
                    % and subject indices
                    Outer_col(idx) = o;
                    MC_col(idx) = mc;
                    N_col(idx) = PopSplits.SampleSizes(n);
                    Subject_col(idx) = Subjects;

                    % Information about construction scheme
                    Topology_col(idx) = string(topology);
                    Weights_col(idx) = string(weights);

                    % Computed central and tail statistics
                    Central_col(idx) = double(reshape(Central(mc,n,:),[],1));
                    Tail_col(idx) = double(reshape(Tail(mc,n,:),[],1));

                    row = row + nr;
                end
            end
        end
    end
end

assert(row == nRows, ...
    'Unexpected number of rows in RepresentativenessStats_AcrossHarmonics!');

% Builds the table per se
RepresentativenessStats_AcrossHarmonics = table( ...
    Outer_col, ...
    MC_col, ...
    N_col, ...
    Subject_col, ...
    Topology_col, ...
    Weights_col, ...
    Central_col, ...
    Tail_col, ...
    'VariableNames',{'Outer','MC','N','Subject','Topology','Weights', ...
    RepCfg.Central.VariableName,RepCfg.Tail.VariableName});

% Makes outer split and Monte Carlo realization categorical (for later use
% in mixed models)
RepresentativenessStats_AcrossHarmonics.Outer = ...
    categorical(RepresentativenessStats_AcrossHarmonics.Outer);

RepresentativenessStats_AcrossHarmonics.MC = ...
    categorical(RepresentativenessStats_AcrossHarmonics.MC);


% RLong concatenates the central and tail summaries into long format; the
% Metric field specifies which summary is stored in Value
RLong = stack(RepresentativenessStats_AcrossHarmonics,RepCfg.MetricNames, ...
    'NewDataVariableName','Value','IndexVariableName','Metric');



%% 3. Global Representativeness means across subjects and realizations
% Computes the mean Representativeness value across subjects, outer splits
% and Monte Carlo realizations for each sample size, construction scheme
% and metric. These values are used for the cross-pillar convergence
% comparison (Figure 5C)
disp('3. Global Representativeness means across subjects and realizations...');

% These are what we will not average over
GroupVars_Global = {'N','Topology','Weights','Metric'};

% Finds all the individual observations from RLong that belong to the same
% group
[G,TGroups] = findgroups(RLong(:,GroupVars_Global));

% Averages within each group
MeanVal = splitapply(@mean,RLong.Value,G);

% Creates the resulting table
RepresentativenessGlobalMeans = [TGroups,table(MeanVal,'VariableNames',{'Mean'})];



%% 4. Generation of cross-subject variability tables
% (a) E_{o,mc}[SD_s(X)]
% Typical across-subject dispersion observed within individual population 
% realizations
% (b) SD_s[E_{o,mc}(X)]
% Across-subject dispersion after averaging across population realizations, 
% i.e., dispersion of subjects' characteristic values
%
% The strategy is to extract from the exhaustive RLong table what we need:
% either (a) all values across subjects for a given {o,mc} realization,
% sample size and construction scheme, to first compute standard deviation
% across subjects and then average across realizations; or (b) all
% realizations available for each subject, to first average across
% realizations and then compute standard deviation across subjects
disp('4. Generation of cross-subject variability tables...');

% The lines below focus on (a)

% Groups corresponding to individual population realizations
GroupVars_Realization = {'Outer','MC','N','Topology','Weights','Metric'};
[G,TRealizationGroups] = findgroups(RLong(:,GroupVars_Realization));

% Computes SD across subjects separately within each realization
SD_AcrossSubjects = splitapply(@std,RLong.Value,G);

% Stores realization-level SD values temporarily in a table
RealizationSD = [TRealizationGroups,table(SD_AcrossSubjects)];

% Groups across population realizations at each sample size
% The idea here is to be able to extract SD values across realizations, so
% that the expectation can be quantified
GroupVars_N = {'N','Topology','Weights','Metric'};
[G,TGroups] = findgroups(RealizationSD(:,GroupVars_N));

% Characterizes the distribution of realization-wise SD values
% xSD contains values across {o,mc}, N, Topology/Weights and Metric
xSD = RealizationSD.SD_AcrossSubjects;

% Mean, 2.5/97.5th percentiles are applied across realizations; they are
% shown in Figure 3B (darker green shade)
% The resulting vectors contain one value per N, construction scheme and
% metric
E_Realization_SD_Subject = splitapply(@mean,xSD,G);
P025_Realization_SD_Subject = splitapply(@(z) prctile(z,2.5),xSD,G);
P975_Realization_SD_Subject = splitapply(@(z) prctile(z,97.5),xSD,G);

% Number of realizations present in each case; we expect n_o * n_MC every
% time
NumRealizations = splitapply(@numel,xSD,G);

% The lines below focus on (b)

% Groups all observations belonging to the same subject at each sample size
% In other words, all {o,mc} realizations for a given
% subject/N/construction scheme are grouped together
GroupVars_Subject = {'Subject','N','Topology','Weights','Metric'};
[G,TSubjectGroups] = findgroups(RLong(:,GroupVars_Subject));

% Computes one characteristic value for each subject by averaging across
% all population realizations in which that subject was left out
MeanAcrossRealizations = splitapply(@mean,RLong.Value,G);

% Here, the number of observations may vary, as a subject can be left out
% in more than one outer split; for each such split, there are NumMC_Run
% values (one per Monte Carlo realization)
NumObservations = splitapply(@numel,RLong.Value,G);

% Contains the mean across realizations, and number of observations,
% for each subject at every N, construction scheme, for each metric of
% interest
% This Table is used to generate the stability of subject ordering (Figure
% 3B, right panel)
RepresentativenessSubjectMeans = [TSubjectGroups,table(NumObservations,...
        MeanAcrossRealizations)];

% Computes across-subject statistics after averaging across realizations
[G,TGroups2] = findgroups(RepresentativenessSubjectMeans(:,GroupVars_N));

% Both aggregation routes must produce the same groups
assert(isequal(TGroups,TGroups2), ...
    'Inconsistent groups across Representativeness variability analyses!');

% Samples each subject's mean across realizations
xSubjectMean = RepresentativenessSubjectMeans.MeanAcrossRealizations;

% Computes mean and standard deviation across subjects for all categories
% (N, construction scheme, metric)
Mean_Subject_E_Realization = splitapply(@mean,xSubjectMean,G);
SD_Subject_E_Realization = splitapply(@std,xSubjectMean,G);

% Number of unique subjects contributing to each estimate
NumUniqueSubjects = splitapply(@numel,xSubjectMean,G);

% We now build the table containing all the relevant data
RepresentativenessSubjectVariability = [TGroups,table(...
    NumRealizations,NumUniqueSubjects, ...
    E_Realization_SD_Subject,P025_Realization_SD_Subject,P975_Realization_SD_Subject, ...
    Mean_Subject_E_Realization,SD_Subject_E_Realization)];

% Descriptive ratio between inter-individual variability retained after
% averaging across population realizations and typical inter-individual
% variability observed within individual realizations
% These values are presented in Supplementary Table 7
%
% NOTE: this is NOT a bounded variance decomposition and may therefore
% exceed one
RepresentativenessSubjectVariability.StabilityRatio = ...
    RepresentativenessSubjectVariability.SD_Subject_E_Realization./...
    RepresentativenessSubjectVariability.E_Realization_SD_Subject;



%% 5. Generation of harmonic-resolved Representativeness table
% Computes harmonic-resolved Representativeness statistics for every
% sample size and connectome construction scheme.
%
% Realization-level summaries across subjects:
%   E_{o,mc}[Mean_s(X)]
%   E_{o,mc}[SD_s(X)]
%
% Subject-level summaries after averaging across realizations:
%   Mean_s[E_{o,mc}(X)]
%   SD_s[E_{o,mc}(X)]
%
% Their ratio is also computed as:
%   SD_s[E_{o,mc}(X)] / E_{o,mc}[SD_s(X)]
%
% One row is generated per N x Harmonic x Topology x Weights.
disp('5. Generation of harmonic-resolved Representativeness table...');

% Number of sample size values
nN = numel(PopSplits.SampleSizes);

% Number of regions/harmonics
nH = Cfg.Atlas.NumRegionsAnalysis;

% Number of subjects
nSubjects = size(Data.Adj,2);

% Number of realizations
nRealizations = NumOuter_Run * NumMC_Run;

% Initializes the table
RepresentativenessHarmonicDescriptives = table;

% Same process for all connectome construction schemes
for it = 1:numel(Topology)

    topology = Topology{it};

    for iw = 1:numel(Weights)

        weights = Weights{iw};

        % (a) Accumulation across population realizations

        % Realization-level mean and SD across subjects
        % Dimensions: realization x N x harmonic
        RealizationMean = nan(nRealizations,nN,nH);
        RealizationSD = nan(nRealizations,nN,nH);

        % Subject-level accumulation across realizations
        % Dimensions: subject x N x harmonic
        SubjectSum = zeros(nSubjects,nN,nH);
        SubjectCount = zeros(nSubjects,nN,nH);

        realization = 0;

        % Accumulation across all outer splits
        for o = 1:NumOuter_Run

            % Held-out subjects for the current outer split
            Subjects = RepMeasSubjects{o};

            % Dimensions: MC x N x subject x harmonic
            X = RepMeas{o}.(topology).(weights);
            
            assert(size(X,1) == NumMC_Run, ...
                'Unexpected number of Monte Carlo realizations in RepMeas!');
            
            % Mean and standard deviation across subjects
            % NaN values are omitted so that across-subject statistics can still be
            % estimated for required-window metrics when the criterion was not reached
            % for some subjects
            MeanSub = mean(X,3,'omitnan');
            SDSub = std(X,[],3,'omitnan');

            % Removes singleton dimension
            % Size is now MC x N x Harmonic
            MeanSub = reshape(MeanSub,NumMC_Run,nN,nH);
            SDSub = reshape(SDSub,NumMC_Run,nN,nH);

            % Indices of the realizations
            idxReal = realization + (1:NumMC_Run);

            % Stores Mean_s(X) and SD_s(X) for all {o,mc} realizations
            % (dimension 1), sample sizes and harmonics
            RealizationMean(idxReal,:,:) = MeanSub;
            RealizationSD(idxReal,:,:) = SDSub;

            % Updates the realization counter (+ NumMC_Run every outer
            % split)
            realization = realization + NumMC_Run;

            % Subject-level accumulation across realizations

            % A subject may be left out in more than one outer split.
            % Values are therefore accumulated according to actual subject
            % identity across all available {o,mc} realizations
            for s = 1:numel(Subjects)

                % Index of currently considered subject
                idSubject = Subjects(s);

                % Dimensions: MC x N x harmonic after singleton removal
                Xs = squeeze(X(:,:,s,:));
                Xs = reshape(Xs,NumMC_Run,nN,nH);

                % Sum and count valid observations across MC realizations
                XsSum = sum(Xs,1,'omitnan');
                XsCount = sum(~isnan(Xs),1);

                % Sum and count across {o,mc} are updated
                SubjectSum(idSubject,:,:) = ...
                    SubjectSum(idSubject,:,:) + XsSum;

                SubjectCount(idSubject,:,:) = ...
                    SubjectCount(idSubject,:,:) + XsCount;
            end
        end

        % Checks that all expected population realizations were included
        assert(realization == nRealizations, ...
            'Unexpected number of Representativeness realizations!');

        % (b) Realization-level quantities

        % Mean across realizations of the mean across subjects
        E_Realization_Mean_Subject = squeeze(mean(RealizationMean,1,'omitnan'));

        % Mean across realizations of the SD across subjects
        E_Realization_SD_Subject = squeeze(mean(RealizationSD,1,'omitnan'));

        % Number of realizations contributing a valid across-subject mean for each
        % N x harmonic. This can vary for required-window metrics if no subject
        % yields a valid value within a given realization
        NumRealizations = squeeze(sum(~isnan(RealizationMean),1));

        % (c) Subject-level quantities after averaging across realizations

        % Characteristic value of each subject after averaging across all
        % population realizations in which that subject was left out
        SubjectMean = SubjectSum./SubjectCount;
        SubjectMean(SubjectCount == 0) = NaN;

        % Mean and SD across subjects
        Mean_Subject_E_Realization = ...
            squeeze(mean(SubjectMean,1,'omitnan'));

        SD_Subject_E_Realization = ...
            squeeze(std(SubjectMean,[],1,'omitnan'));

        % Number of unique subjects contributing to each estimate
        NumUniqueSubjects = squeeze(sum(~isnan(SubjectMean),1));


        % (d) Stability ratio

        % Descriptive ratio between inter-individual variability retained
        % after averaging across population realizations and typical
        % inter-individual variability within individual realizations.
        %
        % NOTE: this is NOT a bounded variance decomposition and may
        % therefore exceed one.
        StabilityRatio = SD_Subject_E_Realization./E_Realization_SD_Subject;
        
        % (e) Conversion to long-format table

        % Creates an N x harmonic grid
        [N_idx,Harmonic_idx] = ndgrid(1:nN,1:nH);

        % After aggregation across realizations and subjects, sample size and
        % harmonic are the remaining dimensions within each construction scheme
        nRows = nN*nH;

        % Values for sample sizes
        % We vectorize the created grid, so we obtain a one-dimensional
        % arrangement of values
        N_col = PopSplits.SampleSizes(N_idx(:));
        N_col = N_col(:);

        % Values for harmonic indices
        Harmonic_col = Harmonic_idx(:);

        % Values for connectome construction scheme
        Topology_col = repmat(string(topology),nRows,1);
        Weights_col = repmat(string(weights),nRows,1);

        NewTable = table(N_col,Harmonic_col,Topology_col,Weights_col, ...
            NumRealizations(:),NumUniqueSubjects(:), ...
            E_Realization_Mean_Subject(:),E_Realization_SD_Subject(:), ...
            Mean_Subject_E_Realization(:),SD_Subject_E_Realization(:), ...
            StabilityRatio(:), ...
            'VariableNames',{'N','Harmonic','Topology','Weights', ...
            'NumRealizations','NumUniqueSubjects', ...
            'E_Realization_Mean_Subject','E_Realization_SD_Subject', ...
            'Mean_Subject_E_Realization','SD_Subject_E_Realization', ...
            'StabilityRatio'});

        RepresentativenessHarmonicDescriptives = ...
            [RepresentativenessHarmonicDescriptives; NewTable];
    end
end



%% 6. Subject-ranking stability across population sample sizes
%
% Quantifies the stability of between-subject ordering in
% Representativeness as population sample size increases.
%
% For each population sample size, subject-level mean Representativeness
% is compared with the corresponding values at N = 700 using Spearman
% correlation.
%
% The analysis uses median matched harmonic similarity for the reference
% AllEdges/AllValues construction
disp('6. Subject-ranking stability across population sample sizes...');

if strcmp(RepCfg.MetricBaseName,'HarmonicSimilarityMatch')

    Topology = "AllEdges";
    Weights = "AllValues";
    Metric = "Median_HarmonicSimilarityMatch";
    
    % Selects the relevant subject-level means
    T = RepresentativenessSubjectMeans( ...
        RepresentativenessSubjectMeans.Topology == Topology & ...
        RepresentativenessSubjectMeans.Weights == Weights & ...
        RepresentativenessSubjectMeans.Metric == Metric,:);
    
    % Population sample sizes
    NSamples = unique(T.N);
    NReference = max(NSamples);
    
    % Reference subject-level values at N = 700
    TRef = T(T.N == NReference,{'Subject','MeanAcrossRealizations'});
    TRef.Properties.VariableNames{'MeanAcrossRealizations'} = 'MeanReference';
    
    % Initializes output
    Rho_N700 = nan(numel(NSamples),1);
    NumSubjects = nan(numel(NSamples),1);
    
    for ns = 1:numel(NSamples)
    
        % Subject-level values at current sample size
        TN = T(T.N == NSamples(ns),{'Subject','MeanAcrossRealizations'});
        TN.Properties.VariableNames{'MeanAcrossRealizations'} = 'MeanCurrent';
    
        % Exact subject matching with the N = 700 reference
        TPair = innerjoin( ...
            TN,TRef, ...
            'Keys','Subject');
    
        % Removes any incomplete pairs
        idx = ...
            isfinite(TPair.MeanCurrent) & ...
            isfinite(TPair.MeanReference);
    
        % Number of subjects entering the correlation
        NumSubjects(ns) = sum(idx);
    
        % Spearman correlation with subject ordering at N = 700
        Rho_N700(ns) = corr( ...
            TPair.MeanCurrent(idx), ...
            TPair.MeanReference(idx), ...
            'Type','Spearman');
    
    end
    
    % Final table
    RepresentativenessSubjectRanking = table( ...
        NSamples, ...
        Rho_N700, ...
        NumSubjects, ...
        'VariableNames',{ ...
            'N', ...
            'Rho_N700', ...
            'NumSubjects'});
    
    assert(RepresentativenessSubjectRanking.N(end) == 700);
    
    assert(abs(RepresentativenessSubjectRanking.Rho_N700(end) - 1) < 1e-12);
end



%% 7. Generation of Linear Mixed Models tables
% Tests changes in Representativeness between predefined sample size
% contrasts for the reference connectome construction scheme.
%
% Delta values are direction-adjusted such that positive values always
% indicate an improvement in Representativeness with increasing sample
% size.
%
% Linear mixed models account for repeated observations across outer
% splits, Monte Carlo realizations and subjects:
%
%   Delta ~ 1 + (1|Outer) + (1|OuterMC) + (1|Subject)
%
% Tangibility is assessed from the model-estimated improvement using
% predefined metric-specific thresholds.
disp('7. Generation of Linear Mixed Models tables...');

% Extracts the reference population-harmonic construction scheme
T = RepresentativenessStats_AcrossHarmonics(...
    RepresentativenessStats_AcrossHarmonics.Topology == RepCfg.ReferenceTopology & ...
    RepresentativenessStats_AcrossHarmonics.Weights == RepCfg.ReferenceWeights,:);

% Sample size contrasts, metrics and corresponding improvement directions
N_Contrasts = RepCfg.N_Contrasts;
Metrics = RepCfg.MetricNames;
ImprovementDirection = RepCfg.ImprovementDirection;

% Identifies the statistic for which a tangibility criterion was
% prespecified, when applicable
if strcmp(RepCfg.TangibilityStatistic,'Central')
    TangibilityMetric = RepCfg.Central.VariableName;
elseif strcmp(RepCfg.TangibilityStatistic,'Tail')
    TangibilityMetric = RepCfg.Tail.VariableName;
else
    TangibilityMetric = "";
end

% Will contain final model outputs
RepresentativenessContrasts = table;

% Same operations for each metric
for id_metric = 1:numel(Metrics)

    metric = Metrics{id_metric};
    direction = ImprovementDirection(id_metric);

    % Same analyses across all contrasts of interest
    for id_contrast = 1:size(N_Contrasts,1)

        % Actual sample sizes at hand
        N1 = N_Contrasts(id_contrast,1);
        N2 = N_Contrasts(id_contrast,2);

        % Extracts the metrics data associated to both sample sizes
        T1 = T(T.N == N1,:);
        T2 = T(T.N == N2,:);

        % Ensures identical ordering across paired trajectories
        T1 = sortrows(T1,{'Outer','MC','Subject'});
        T2 = sortrows(T2,{'Outer','MC','Subject'});

        % Checks that observations are correctly paired
        assert(height(T1) == height(T2), ...
            'Different numbers of observations across paired sample sizes!');
        assert(isequal(T1.Outer,T2.Outer), ...
            'Outer splits do not match across paired sample sizes!');
        assert(isequal(T1.MC,T2.MC), ...
            'Monte Carlo realizations do not match across paired sample sizes!');
        assert(isequal(T1.Subject,T2.Subject), ...
            'Subjects do not match across paired sample sizes!');

        % Computes direction-adjusted changes so that positive Delta always
        % indicates improvement in Representativeness
        if direction == 1
            Delta = T2.(metric) - T1.(metric);
        else
            Delta = T1.(metric) - T2.(metric);
        end

        % Builds the statistical table for the current contrast
        TDelta = table(categorical(T1.Outer),categorical(T1.MC), ...
            categorical(string(T1.Outer) + "_" + string(T1.MC)), ...
            categorical(T1.Subject), ...
            Delta, ...
            'VariableNames',{'Outer','MC','OuterMC','Subject','Delta'});

        % Fits the linear mixed model:
        % Delta_oms = beta_0 + u_o + v_om + w_s + epsilon_oms
        mdl = fitlme(TDelta, ...
            'Delta ~ 1 + (1|Outer) + (1|OuterMC) + (1|Subject)');

        % Extracts fixed-effect statistics for the intercept, corresponding
        % to the model-estimated direction-adjusted change
        C = mdl.Coefficients;

        Estimate = C.Estimate(1);
        SE = C.SE(1);
        DF = C.DF(1);
        tStat = C.tStat(1);
        pValue = C.pValue(1);
        CI_Lower = C.Lower(1);
        CI_Upper = C.Upper(1);

        % Extracts random-effect and residual variance components
        [Psi,MSE] = covarianceParameters(mdl);

        Var_Outer = Psi{1};
        Var_OuterMC = Psi{2};
        Var_Subject = Psi{3};
        Var_Residual = MSE;

        % Empirical Delta summaries reported as
        % median [2.5th percentile, 97.5th percentile]
        Median_Delta = median(Delta);
        P025_Delta = prctile(Delta,2.5);
        P975_Delta = prctile(Delta,97.5);

        % Determines whether the model-estimated improvement reaches the
        % predefined metric-specific threshold for a tangible change
        % Tangibility is assessed only for the prespecified statistic
        if strcmp(metric,TangibilityMetric)
        
            TangibilityThreshold = RepCfg.TangibilityThreshold;
            Tangible = Estimate >= TangibilityThreshold;
        
        else
        
            TangibilityThreshold = NaN;
            Tangible = NaN;
        
        end

        % Contrast identifier
        Contrast = string(N1) + "_to_" + string(N2);

        % Adds the current contrast to the final summary table
        NewRow = table(string(metric),N1,N2,Contrast, ...
            Median_Delta,P025_Delta,P975_Delta, ...
            Estimate,CI_Lower,CI_Upper,SE,DF,tStat,pValue, ...
            Var_Outer,Var_OuterMC,Var_Subject,Var_Residual, ...
            TangibilityThreshold,Tangible, ...
            'VariableNames',{'Metric','N1','N2','Contrast', ...
            'Median_Delta','P025_Delta','P975_Delta', ...
            'Estimate','CI_Lower','CI_Upper','SE','DF','tStat','pValue', ...
            'Var_Outer','Var_OuterMC','Var_Subject','Var_Residual', ...
            'TangibilityThreshold','Tangible'});

        RepresentativenessContrasts = ...
            [RepresentativenessContrasts; NewRow];
    end
end



%% 8. Saving all Representativeness results
% Saves all outputs associated with the Representativeness metric
% considered in the current run
disp('8. Saving all Representativeness results...');

OutputFile = sprintf('Representativeness_%s_Results.mat', ...
    RepCfg.MetricBaseName);

save(fullfile(ResultsDir,OutputFile), ...
    'RepCfg', ...
    'RepMeas', ...
    'RepMeasSubjects', ...
    'RepresentativenessStats_AcrossHarmonics', ...
    'RepresentativenessGlobalMeans', ...
    'RepresentativenessSubjectMeans', ...
    'RepresentativenessSubjectVariability', ...
    'RepresentativenessHarmonicDescriptives', ...
    'RepresentativenessContrasts', ...
    '-v7.3');

% Subject-ranking stability is computed only for matched harmonic similarity
if strcmp(RepCfg.MetricBaseName,'HarmonicSimilarityMatch')
    save(fullfile(ResultsDir,OutputFile), ...
        'RepresentativenessSubjectRanking', ...
        '-append');
end