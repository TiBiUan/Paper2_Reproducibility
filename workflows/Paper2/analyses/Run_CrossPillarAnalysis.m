%% 1. Robustness versus Representativeness
%
% Examines whether population representations that are more robust across
% independently sampled populations are also more representative of
% subject-specific harmonic organization.
%
% Analyses use the reference AllEdges/AllValues construction.
%
% Robustness is quantified using median matched harmonic similarity between
% independently sampled population representations.
%
% Representativeness is quantified using median matched harmonic similarity
% between population and subject-specific representations. Subject-level
% values are averaged within each {Outer,MC,N} realization to provide one
% Representativeness value corresponding to each Robustness realization.
%
% Two complementary associations are assessed:
%
%   1. Spearman correlation between Robustness and Representativeness
%      separately at each population sample size.
%
%   2. A pooled sample size-adjusted association obtained by rank-
%      transforming Robustness and Representativeness separately within
%      each sample size and correlating the resulting normalized ranks.

% Reference construction is considered for the analyses below

load(fullfile(ResultsDir,'Robustness_Results.mat'),'RobustnessStats');
load(fullfile(ResultsDir,...
    'Representativeness_HarmonicSimilarityMatch_Results.mat'),...
    'RepresentativenessStats_AcrossHarmonics');

% Robustness
T_Rob = RobustnessStats(RobustnessStats.Topology == "AllEdges" & ...
    RobustnessStats.Weights == "AllValues",:);

% Representativeness
T_Rep = RepresentativenessStats_AcrossHarmonics(...
    RepresentativenessStats_AcrossHarmonics.Topology == "AllEdges" & ...
    RepresentativenessStats_AcrossHarmonics.Weights == "AllValues",:);

% Average Representativeness across subjects (and SD too) is computed for 
% each outer split, Monte Carlo realization and sample size
GroupVars = {'Outer','MC','N'};
[G,TGroups] = findgroups(T_Rep(:,GroupVars));

MeanRepresentativeness = splitapply(@mean,...
    T_Rep.Median_HarmonicSimilarityMatch,G);

SDRepresentativeness = splitapply(@std,...
    T_Rep.Median_HarmonicSimilarityMatch,G);

% Puts the results into a table
Rep_AvgSubj = [TGroups,table(MeanRepresentativeness,SDRepresentativeness)];

% Matches Robustness and Representativeness realizations

% Quantities required from Robustness
Rob_Link = T_Rob(:,{'Outer','MC','N','Median_HarmSim_Index'});
Rob_Link.Properties.VariableNames{'Median_HarmSim_Index'} = 'Robustness';

% Quantities required from Representativeness
Rep_Link = Rep_AvgSubj(:,{'Outer','MC','N','MeanRepresentativeness','SDRepresentativeness'});

% Explicitly matches identical {Outer,MC,N} realizations
RobustnessRepresentativenessStats = innerjoin(Rob_Link,Rep_Link, ...
    'Keys',{'Outer','MC','N'});

% Ensures that all realizations were matched uniquely
assert(height(RobustnessRepresentativenessStats) == height(Rob_Link));

% Sample size-specific associations

% All sample sizes of interest
NSamples = unique(RobustnessRepresentativenessStats.N);

% Will contain sample size-wise associations
SpearmanRho = nan(numel(NSamples),1);

% Number of realizations over which we compute the association for each
% sample size
NumRealizations = nan(numel(NSamples),1);

% Same process for each sample size
for ns = 1:numel(NSamples)

    % Current sample size
    N = NSamples(ns);

    % Corresponding realizations
    idx = RobustnessRepresentativenessStats.N == N;

    % Number of matched realizations
    NumRealizations(ns) = sum(idx);

    % Spearman correlation
    SpearmanRho(ns) = corr(RobustnessRepresentativenessStats.Robustness(idx), ...
        RobustnessRepresentativenessStats.MeanRepresentativeness(idx), ...
        'Type','Spearman','Rows','complete');
end

% Summary across sample sizes
RobustnessRepresentativenessByN = table(NSamples,NumRealizations,SpearmanRho, ...
    'VariableNames',{'N','NumRealizations','SpearmanRho'});

% Pooled sample size-adjusted association

% Converts Robustness and Representativeness to normalized ranks separately
% within each population sample size.
RankRobustness = nan(height(RobustnessRepresentativenessStats),1);
RankRepresentativeness = nan(height(RobustnessRepresentativenessStats),1);

for ns = 1:numel(NSamples)

    % Current sample size
    idx = RobustnessRepresentativenessStats.N == NSamples(ns);

    % Number of realizations
    n = sum(idx);

    % Normalized within-N ranks
    RankRobustness(idx) = ...
        (tiedrank(RobustnessRepresentativenessStats.Robustness(idx))-0.5)/n;

    RankRepresentativeness(idx) = ...
        (tiedrank( ...
        RobustnessRepresentativenessStats.MeanRepresentativeness(idx))-0.5)/n;
end

% Stores normalized ranks with matched observations
RobustnessRepresentativenessStats.RankRobustness = RankRobustness;
RobustnessRepresentativenessStats.RankRepresentativeness = RankRepresentativeness;

% Pooled association after adjustment for sample size
R_PooledWithinN = corr(RankRobustness,RankRepresentativeness, ...
    'Type','Pearson','Rows','complete');



%% 2. Representativeness versus Practical Utility
%
% Examines whether subjects that are more representative of the population
% harmonic basis also show smaller errors when population-derived rather
% than subject-specific harmonics are used for functional analyses.
%
% Representativeness is quantified using median matched harmonic similarity.
%
% Practical Utility is quantified using:
%
%   StructuralMismatch    Structural mismatch
%   Error_Functional      Functional error
%   Rel_DeltaC            Relative second-order moment error
%   Error_FC              FC error
%
% Subject-level Representativeness and Practical Utility values are matched
% within each {Outer,MC,N} realization.
%
% Spearman correlations are computed separately for each realization,
% Practical Utility metric and graph-frequency range.
%
% Correlations are Fisher Z-transformed for statistical analyses.

% Reference construction is considered for the analyses below

load(fullfile(ResultsDir,...
    'Representativeness_HarmonicSimilarityMatch_Results.mat'),...
    'RepresentativenessStats_AcrossHarmonics');

load(fullfile(ResultsDir,'PracticalUtility_Part2_Results.mat'),'PracticalUtilitySSStats');

% Representativeness
T_Rep = RepresentativenessStats_AcrossHarmonics( ...
    RepresentativenessStats_AcrossHarmonics.Topology == "AllEdges" & ...
    RepresentativenessStats_AcrossHarmonics.Weights == "AllValues",:);

T_Rep = T_Rep(:,{'Outer','MC','N','Subject','Median_HarmonicSimilarityMatch'});

% Practical Utility
T_PU = PracticalUtilitySSStats(:,{'Outer','MC','N','Subject','Filter', ...
    'StructuralMismatch','Error_Functional','Rel_DeltaC','Error_FC'});

% Practical Utility metrics of interest
PUMetrics = {'StructuralMismatch','Error_Functional','Rel_DeltaC','Error_FC'};

% Filter sub-types
Filters = {'Low','High'};

% Matches subject-level Representativeness and Practical Utility

% Explicitly matches identical {Outer,MC,N,Subject} observations
RepresentativenessPracticalUtilityStats = innerjoin(T_Rep,T_PU, ...
    'Keys',{'Outer','MC','N','Subject'});

% Ensures that every Practical Utility observation was matched
assert(height(RepresentativenessPracticalUtilityStats) == height(T_PU));

% Computes correlations within each realization

% Will apply the same for all {o,mc} realizations, sample sizes and filter
% sub-types (i.e., an operation across subjects)
GroupVars = {'Outer','MC','N','Filter'};
[G,TGroups] = findgroups(RepresentativenessPracticalUtilityStats(:,GroupVars));

% Number of realization × filter combinations
nGroups = height(TGroups);

% Number of Practical Utility metrics
nMetrics = numel(PUMetrics);

% Initializes long-form results
Outer = repmat(TGroups.Outer,nMetrics,1);
MC = repmat(TGroups.MC,nMetrics,1);
N = repmat(TGroups.N,nMetrics,1);
Filter = repmat(TGroups.Filter,nMetrics,1);

Metric = strings(nGroups*nMetrics,1);
SpearmanRho = nan(nGroups*nMetrics,1);
FisherZ = nan(nGroups*nMetrics,1);
NumSubjects = nan(nGroups*nMetrics,1);

% Avoids numerical problems if rho happens to be exactly +/-1
cliprho = @(r) min(max(r,-1+eps),1-eps);

% Same process for each Practical Utility metric
for m = 1:nMetrics

    % Current metric
    metric = PUMetrics{m};

    % Rows corresponding to this metric in the long-form output
    rows = (m-1)*nGroups + (1:nGroups);

    % Metric name
    Metric(rows) = string(metric);

    % Spearman correlation within each {Outer,MC,N,Filter}
    SpearmanRho(rows) = splitapply(@(x,y) corr(x,y,'Type','Spearman','Rows','complete'), ...
        RepresentativenessPracticalUtilityStats.Median_HarmonicSimilarityMatch, ...
        RepresentativenessPracticalUtilityStats.(metric),G);

    % Number of subjects contributing to each correlation
    NumSubjects(rows) = splitapply(@(x,y) sum(~isnan(x) & ~isnan(y)), ...
        RepresentativenessPracticalUtilityStats.Median_HarmonicSimilarityMatch, ...
        RepresentativenessPracticalUtilityStats.(metric),G);

    % Fisher Z transformation
    FisherZ(rows) = atanh(cliprho(SpearmanRho(rows)));
end

% Correlation table
RepresentativenessPracticalUtilityCorrelations = table(Outer,MC,N,Filter,...
    Metric,NumSubjects, ...
    SpearmanRho,FisherZ);

% Low-versus-High frequency comparison

% Separates Low and High correlations
T_Low = RepresentativenessPracticalUtilityCorrelations( ...
    RepresentativenessPracticalUtilityCorrelations.Filter == "Low",:);

T_High = RepresentativenessPracticalUtilityCorrelations( ...
    RepresentativenessPracticalUtilityCorrelations.Filter == "High",:);

% Quantities required for matching
T_Low = T_Low(:,{ ...
    'Outer','MC','N','Metric','SpearmanRho','FisherZ'});

T_High = T_High(:,{ ...
    'Outer','MC','N','Metric','SpearmanRho','FisherZ'});

% Makes frequency explicit in variable names
T_Low.Properties.VariableNames{'SpearmanRho'} = 'Rho_Low';
T_Low.Properties.VariableNames{'FisherZ'} = 'Z_Low';

T_High.Properties.VariableNames{'SpearmanRho'} = 'Rho_High';
T_High.Properties.VariableNames{'FisherZ'} = 'Z_High';

% Explicitly matches Low and High observations
RepresentativenessPracticalUtilityFrequencyComparison = innerjoin(T_Low,T_High, ...
    'Keys',{'Outer','MC','N','Metric'});

% Ensures complete one-to-one matching
assert(height(RepresentativenessPracticalUtilityFrequencyComparison) == height(T_Low));
assert(height(RepresentativenessPracticalUtilityFrequencyComparison) == height(T_High));

% Difference between High- and Low-frequency correlations
RepresentativenessPracticalUtilityFrequencyComparison.DeltaZ = ...
    RepresentativenessPracticalUtilityFrequencyComparison.Z_High - ...
    RepresentativenessPracticalUtilityFrequencyComparison.Z_Low;

% Statistical analyses

% Initializes result tables
RepresentativenessPracticalUtilityModels = table;
RepresentativenessPracticalUtilityFrequencyModels = table;

% Correlation models

% Same process for each Practical Utility metric and filter
for m = 1:nMetrics

    metric = PUMetrics{m};

    for id_filter = 1:numel(Filters)

        filter = Filters{id_filter};

        % Current observations
        T = RepresentativenessPracticalUtilityCorrelations( ...
            RepresentativenessPracticalUtilityCorrelations.Metric == metric & ...
            RepresentativenessPracticalUtilityCorrelations.Filter == filter,:);

        % Random-effect identifiers
        T.Outer = categorical(T.Outer);

        T.OuterMC = categorical(string(T.Outer) + "_" + string(T.MC));

        % Centered log sample size
        T.logN = log(T.N);
        T.logNc = T.logN - mean(T.logN);

        % Mixed-effects model
        mdl = fitlme(T, ...
            'FisherZ ~ 1 + logNc + (1|Outer) + (1|OuterMC)');

        % Fixed effects
        C = mdl.Coefficients;

        % Intercept: average Fisher-Z association at centered log(N)
        row = table(string(metric),string(filter), ...
            "Intercept", ...
            C.Estimate(1),C.SE(1),C.DF(1),C.tStat(1),C.pValue(1),C.Lower(1),C.Upper(1), ...
            'VariableNames',{'Metric','Filter',...
            'Effect','Estimate','SE','DF','tStat','pValue','CI_Lower','CI_Upper'});

        RepresentativenessPracticalUtilityModels = ...
            [RepresentativenessPracticalUtilityModels;row];

        % log(N) effect
        row = table(string(metric),string(filter), ...
            "logN", ...
            C.Estimate(2),C.SE(2),C.DF(2),C.tStat(2),C.pValue(2),C.Lower(2),C.Upper(2), ...
            'VariableNames',{'Metric','Filter',...
            'Effect','Estimate','SE','DF','tStat','pValue','CI_Lower','CI_Upper'});

        RepresentativenessPracticalUtilityModels = ...
            [RepresentativenessPracticalUtilityModels;row];
    end
end

% Low-versus-High frequency models

% Same process for each Practical Utility metric
for m = 1:nMetrics

    metric = PUMetrics{m};

    % Structural mismatch is identical for Low and High filtering.
    % Consequently, the Low-High difference in Fisher-Z correlations is
    % identically zero and no frequency-difference model is fitted
    if ~strcmp(metric,'StructuralMismatch')

        % Current observations
        T = RepresentativenessPracticalUtilityFrequencyComparison( ...
            RepresentativenessPracticalUtilityFrequencyComparison.Metric == metric,:);
    
        % Random-effect identifiers
        T.Outer = categorical(T.Outer);
    
        T.OuterMC = categorical( ...
            string(T.Outer) + "_" + string(T.MC));
    
        % Centered log sample size
        T.logN = log(T.N);
        T.logNc = T.logN - mean(T.logN);
    
        % Mixed-effects model
        mdl = fitlme(T, ...
            'DeltaZ ~ 1 + logNc + (1|Outer) + (1|OuterMC)');
    
        % Fixed effects
        C = mdl.Coefficients;
    
        % Stores both intercept and sample size effect
        for c = 1:height(C)
    
            row = table(string(metric),string(C.Name{c}), ...
                C.Estimate(c),C.SE(c),C.DF(c),C.tStat(c),C.pValue(c),C.Lower(c),C.Upper(c), ...
                'VariableNames',{'Metric','Effect', ...
                'Estimate','SE','DF','tStat','pValue','CI_Lower','CI_Upper'});
    
            RepresentativenessPracticalUtilityFrequencyModels = ...
                [RepresentativenessPracticalUtilityFrequencyModels;row];
        end
    end
end



%% 3. Cross-pillar convergence
%
% Compares the evolution of Robustness, Representativeness and Practical
% Utility with increasing population sample size.
%
% For each quantity, convergence is expressed relative to the total change
% observed between N = 5 and N = 350:
%
%                 Q(N) - Q(5)
%   C(N) = -------------------------
%                Q(350) - Q(5)
%
% Consequently, C(5) = 0 and C(350) = 1 for every quantity. Values can
% exceed the [0,1] interval if the trajectory is non-monotonic.
%
% Analyses use the reference AllEdges/AllValues construction and sample
% sizes up to N = 350, corresponding to the range available for all three
% pillars.

% Quantities included in the comparison

% Robustness
RobustnessMetrics = {'Median_HarmSim_Match', ...
    'P95_ProjErr_ReqWin', ...
    'LaplacianReconstructionError'};

% Representativeness
RepresentativenessMetrics = {'Median_HarmonicSimilarityMatch', ...
    'Median_LeakageEntropy', ...
    'Median_BandwidthIndex', ...
    'P95_SubspaceSimilarityRequiredWindows'};

% Practical Utility
PracticalUtilityMetrics = {'StructuralMismatch', ...
    'Error_Functional', ...
    'Rel_DeltaC', ...
    'Error_FC'};

% Robustness
load(fullfile(ResultsDir,'Robustness_Results.mat'),'RobustnessDescriptives');

% Reference construction and common sample size range
T = RobustnessDescriptives(RobustnessDescriptives.Topology == "AllEdges" & ...
    RobustnessDescriptives.Weights == "AllValues" & ...
    RobustnessDescriptives.N <= 350,:);

CrossPillarConvergence = table;

for m = 1:numel(RobustnessMetrics)

    % Current metric
    metric = RobustnessMetrics{m};

    % Corresponding descriptive values
    T_Metric = T(T.Metric == metric,:);
    T_Metric = sortrows(T_Metric,'N');

    % Validates common convergence range
    assert(T_Metric.N(1) == 5);
    assert(T_Metric.N(end) == 350);

    % Mean quantity at each sample size
    Q = T_Metric.Mean;

    % Normalized observed convergence
    Convergence = (Q-Q(1))./(Q(end)-Q(1));

    % Stores results
    T_Out = table(repmat("Robustness",height(T_Metric),1), ...
        repmat(string(metric),height(T_Metric),1), ...
        repmat("None",height(T_Metric),1), ...
        T_Metric.N,Q,Convergence, ...
        'VariableNames',{'Pillar', ...
        'Metric', ...
        'Filter', ...
        'N','Mean','Convergence'});

    CrossPillarConvergence = ...
        [CrossPillarConvergence;T_Out];
end

% Representativeness
RepFiles = {'Representativeness_HarmonicSimilarityMatch_Results.mat'
    'Representativeness_LeakageEntropy_Results.mat'
    'Representativeness_BandwidthIndex_Results.mat'
    'Representativeness_SubspaceSimilarityRequiredWindows_Results.mat'};

for m = 1:numel(RepFiles)

    % Loads the data for the metric of interest
    RepTmp = load(RepFiles{m},'RepresentativenessGlobalMeans');

    % Gets only the values of interest
    T = RepTmp.RepresentativenessGlobalMeans;
    T = T(T.Topology == "AllEdges" & T.Weights == "AllValues" & T.N <= 350,:);

    % Current metric
    metric = RepresentativenessMetrics{m};

    % Corresponding descriptive values
    T_Metric = T(T.Metric == metric,:);
    T_Metric = sortrows(T_Metric,'N');

    % Validates common convergence range
    assert(T_Metric.N(1) == 5);
    assert(T_Metric.N(end) == 350);

    % Mean quantity at each sample size
    Q = T_Metric.Mean;

    % Normalized observed convergence
    Convergence = (Q-Q(1))./(Q(end)-Q(1));

    % Stores results
    T_Out = table(repmat("Representativeness",height(T_Metric),1), ...
        repmat(string(metric),height(T_Metric),1), ...
        repmat("None",height(T_Metric),1), ...
        T_Metric.N,Q,Convergence, ...
        'VariableNames',{'Pillar', ...
        'Metric', ...
        'Filter', ...
        'N','Mean','Convergence'});

    CrossPillarConvergence = ...
        [CrossPillarConvergence;T_Out];
end

% Practical Utility
load('PracticalUtility_Part2_Results.mat','PracticalUtilitySSDescriptives');

% Common sample size range
T = PracticalUtilitySSDescriptives(PracticalUtilitySSDescriptives.N <= 350,:);

% Filer sub-types
Filters = {'Low','High'};

for m = 1:numel(PracticalUtilityMetrics)

    % Current metric
    metric = PracticalUtilityMetrics{m};

    for id_filter = 1:numel(Filters)

        % Current filter
        filter = Filters{id_filter};

        % Corresponding descriptive values
        T_Metric = T(T.Metric == metric & T.Filter == filter,:);
        T_Metric = sortrows(T_Metric,'N');

        % Validates common convergence range
        assert(T_Metric.N(1) == 5);
        assert(T_Metric.N(end) == 350);

        % Mean quantity at each sample size
        Q = T_Metric.Mean;

        % Normalized observed convergence
        Convergence = (Q-Q(1))./(Q(end)-Q(1));

        % Stores results
        T_Out = table(repmat("PracticalUtility",height(T_Metric),1), ...
            repmat(string(metric),height(T_Metric),1), ...
            repmat(string(filter),height(T_Metric),1), ...
            T_Metric.N,Q,Convergence, ...
            'VariableNames',{ ...
            'Pillar', ...
            'Metric', ...
            'Filter', ...
            'N','Mean','Convergence'});

        CrossPillarConvergence = ...
            [CrossPillarConvergence;T_Out];
    end
end

% Validation of common sample sizes

% Sample sizes expected for the cross-pillar comparison
N_Common = unique(CrossPillarConvergence.N);

% Each trajectory must contain every common sample size exactly once
[G,TrajectoryGroups] = findgroups( ...
    CrossPillarConvergence(:,{'Pillar','Metric','Filter'}));

NumN = splitapply(@numel,CrossPillarConvergence.N,G);

assert(all(NumN == numel(N_Common)), ...
    'Cross-pillar convergence trajectories contain different numbers of sample sizes.');

for g = 1:height(TrajectoryGroups)

    idx = G == g;

    assert(isequal( ...
        sort(CrossPillarConvergence.N(idx)), ...
        sort(N_Common)), ...
        'Cross-pillar convergence trajectories use different sample sizes.');
end

% Every trajectory must start at 0 and end at 1
for g = 1:height(TrajectoryGroups)

    T_Check = CrossPillarConvergence(G == g,:);
    T_Check = sortrows(T_Check,'N');

    assert(abs(T_Check.Convergence(1)) < 1e-12);
    assert(abs(T_Check.Convergence(end)-1) < 1e-12);
end



%% 4. Save results
% Saves all outputs required for downstream analyses and figure generation

save(fullfile(ResultsDir,'CrossPillarAnalysis_Results.mat'), ...
    'RobustnessRepresentativenessStats', ...
    'RobustnessRepresentativenessByN', ...
    'R_PooledWithinN', ...
    'RepresentativenessPracticalUtilityStats', ...
    'RepresentativenessPracticalUtilityCorrelations', ...
    'RepresentativenessPracticalUtilityFrequencyComparison', ...
    'RepresentativenessPracticalUtilityModels', ...
    'RepresentativenessPracticalUtilityFrequencyModels', ...
    'CrossPillarConvergence', ...
    '-v7.3');