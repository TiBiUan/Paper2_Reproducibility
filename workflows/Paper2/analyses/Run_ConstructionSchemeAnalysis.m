%% 1. Robustness across construction schemes
%
% Examines the influence of connectome construction on Robustness.
%
% Two construction factors are considered:
%
%   Topology:
%       AllEdges
%       ConsensusEdges
%
%   Weights:
%       AllValues
%       NZValues
%
% Three complementary analyses are performed:
%
%   1. Construction effects on Robustness at the largest population sample
%      size available for Robustness (N = 350).
%
%   2. Construction effects on Robustness gains across predefined sample
%      size contrasts.
%
%   3. Sensitivity of harmonic-wise Robustness profiles to changes in
%      topology and weight aggregation at N = 350.



%% Endpoint construction effects

load(fullfile(ResultsDir,'Robustness_Results.mat'),'RobustnessStats');
load(fullfile(ResultsDir,'Robustness_Results.mat'),'Robustness');

% Largest population sample size available for Robustness
Nmax_Rob = 350;

% Ensures construction factors are categorical
RobustnessStats.Topology = categorical(RobustnessStats.Topology);
RobustnessStats.Weights = categorical(RobustnessStats.Weights);

% Nested realization identifier
RobustnessStats.OuterMC = categorical( ...
    "O" + string(RobustnessStats.Outer) + ...
    "_MC" + string(RobustnessStats.MC));

% Endpoint observations
T350 = RobustnessStats(RobustnessStats.N == Nmax_Rob,:);

% Robustness outcome
Outcome = 'Median_HarmSim_Match';

% Factorial mixed-effects model
RobustnessConstructionEndpointModel = fitlme(T350, ...
    [Outcome ' ~ Topology*Weights + (1|Outer) + (1|OuterMC)']);

% Simple construction effects
R1 = GetLMEContrast( ...
    RobustnessConstructionEndpointModel,[0 1 0 0]);

R2 = GetLMEContrast( ...
    RobustnessConstructionEndpointModel,[0 1 0 1]);

R3 = GetLMEContrast( ...
    RobustnessConstructionEndpointModel,[0 0 1 0]);

R4 = GetLMEContrast( ...
    RobustnessConstructionEndpointModel,[0 0 1 1]);

% Combines endpoint construction effects
RobustnessConstructionEndpoint = [R1;R2;R3;R4];

RobustnessConstructionEndpoint.Effect = categorical({ ...
    'Topology_AllValues'; ...
    'Topology_NZValues'; ...
    'Weights_AllEdges'; ...
    'Weights_ConsensusEdges'});



%% Construction effects on sample size contrasts

% Predefined sample size contrasts
NContrasts = [ ...
     10   50
     50  100
    100  200
    200  350];

% Will contain matched observations for all contrasts
RobustnessConstructionContrastStats = table();

% Constructs sample size gains
for nc = 1:size(NContrasts,1)

    % Current contrast
    N1 = NContrasts(nc,1);
    N2 = NContrasts(nc,2);

    % Observations at the two sample sizes
    T1 = RobustnessStats( ...
        RobustnessStats.N == N1, ...
        {'Outer','MC','Topology','Weights',Outcome});

    T2 = RobustnessStats( ...
        RobustnessStats.N == N2, ...
        {'Outer','MC','Topology','Weights',Outcome});

    % Renames outcome before matching
    T1.Properties.VariableNames{Outcome} = 'Value_N1';
    T2.Properties.VariableNames{Outcome} = 'Value_N2';

    % Explicitly matches identical realizations and construction schemes
    TD = innerjoin(T1,T2, ...
        'Keys',{'Outer','MC','Topology','Weights'});

    % Ensures complete one-to-one matching
    assert(height(TD) == height(T1));
    assert(height(TD) == height(T2));

    % Larger HWSM indicates greater Robustness, so positive values indicate
    % improvement with increasing population sample size
    TD.Gain = TD.Value_N2 - TD.Value_N1;

    % Stores contrast information
    TD.N1 = repmat(N1,height(TD),1);
    TD.N2 = repmat(N2,height(TD),1);
    TD.Contrast = repmat(nc,height(TD),1);

    % Nested realization identifier
    TD.OuterMC = categorical( ...
        "O" + string(TD.Outer) + ...
        "_MC" + string(TD.MC));

    % Appends observations
    RobustnessConstructionContrastStats = ...
        [RobustnessConstructionContrastStats;TD];
end

% Initializes model and result containers
RobustnessConstructionContrastModels = ...
    cell(size(NContrasts,1),1);

ContrastResults = ...
    cell(size(NContrasts,1),1);

% Fits one factorial model for each sample size contrast
for nc = 1:size(NContrasts,1)

    % Current contrast
    TD = RobustnessConstructionContrastStats( ...
        RobustnessConstructionContrastStats.Contrast == nc,:);

    % Factorial mixed-effects model
    mdl = fitlme(TD, ...
        'Gain ~ Topology*Weights + (1|Outer) + (1|OuterMC)');

    RobustnessConstructionContrastModels{nc} = mdl;

    % Simple construction effects
    R1 = GetLMEContrast(mdl,[0 1 0 0]);
    R2 = GetLMEContrast(mdl,[0 1 0 1]);
    R3 = GetLMEContrast(mdl,[0 0 1 0]);
    R4 = GetLMEContrast(mdl,[0 0 1 1]);

    R = [R1;R2;R3;R4];

    % Labels construction effects
    R.Effect = categorical({ ...
        'Topology_AllValues'; ...
        'Topology_NZValues'; ...
        'Weights_AllEdges'; ...
        'Weights_ConsensusEdges'});

    % Stores contrast information
    R.N1 = repmat(NContrasts(nc,1),4,1);
    R.N2 = repmat(NContrasts(nc,2),4,1);
    R.Contrast = repmat(nc,4,1);

    ContrastResults{nc} = R;
end

% Combines results across sample size contrasts
RobustnessConstructionContrasts = ...
    vertcat(ContrastResults{:});

RobustnessConstructionContrasts = movevars( ...
    RobustnessConstructionContrasts, ...
    {'Contrast','N1','N2','Effect'}, ...
    'Before',1);



%% Spectral-profile sensitivity to construction scheme

% Finds N = 350 in the original Robustness structure
SampleSizes_Rob = Cfg.Population.RobustnessSampleSizes;
iN = find(SampleSizes_Rob == Nmax_Rob,1);

assert(~isempty(iN), ...
    'N = 350 was not found in the Robustness sample sizes!');

% Construction schemes
TopologyNames = {'AllEdges','ConsensusEdges'};
WeightsNames = {'AllValues','NZValues'};

% Will contain one observation per Outer x MC realization
RobustnessProfileSensitivity = table();

for o = 1:numel(Robustness.Outer)

    for mc = 1:numel(Robustness.Outer(o).MonteCarlo)

        % Harmonic-wise Robustness profiles for the four construction
        % schemes
        P_AE_AV = ...
            Robustness.Outer(o).MonteCarlo(mc).SampleSize(iN). ...
            AllEdges.AllValues.HarmonicSimilarity.Match(:);

        P_AE_NZ = ...
            Robustness.Outer(o).MonteCarlo(mc).SampleSize(iN). ...
            AllEdges.NZValues.HarmonicSimilarity.Match(:);

        P_CE_AV = ...
            Robustness.Outer(o).MonteCarlo(mc).SampleSize(iN). ...
            ConsensusEdges.AllValues.HarmonicSimilarity.Match(:);

        P_CE_NZ = ...
            Robustness.Outer(o).MonteCarlo(mc).SampleSize(iN). ...
            ConsensusEdges.NZValues.HarmonicSimilarity.Match(:);

        % Validates harmonic profiles
        assert(numel(P_AE_AV) == Cfg.Atlas.NumRegionsAnalysis);
        assert(numel(P_AE_NZ) == Cfg.Atlas.NumRegionsAnalysis);
        assert(numel(P_CE_AV) == Cfg.Atlas.NumRegionsAnalysis);
        assert(numel(P_CE_NZ) == Cfg.Atlas.NumRegionsAnalysis);

        % Topology sensitivity:
        % changes Topology while keeping Weights fixed
        R_Topology_AllValues = corr( ...
            P_AE_AV,P_CE_AV, ...
            'Type','Pearson', ...
            'Rows','complete');

        R_Topology_NZValues = corr( ...
            P_AE_NZ,P_CE_NZ, ...
            'Type','Pearson', ...
            'Rows','complete');

        % Weight sensitivity:
        % changes Weights while keeping Topology fixed
        R_Weights_AllEdges = corr( ...
            P_AE_AV,P_AE_NZ, ...
            'Type','Pearson', ...
            'Rows','complete');

        R_Weights_ConsensusEdges = corr( ...
            P_CE_AV,P_CE_NZ, ...
            'Type','Pearson', ...
            'Rows','complete');

        % Avoids numerical problems if a correlation equals +/-1
        cliprho = @(r) min(max(r,-1+eps),1-eps);

        % Fisher Z transformation
        Z_Topology_AllValues = ...
            atanh(cliprho(R_Topology_AllValues));

        Z_Topology_NZValues = ...
            atanh(cliprho(R_Topology_NZValues));

        Z_Weights_AllEdges = ...
            atanh(cliprho(R_Weights_AllEdges));

        Z_Weights_ConsensusEdges = ...
            atanh(cliprho(R_Weights_ConsensusEdges));

        % Stores one observation per {Outer,MC}
        tmp = table( ...
            o, ...
            mc, ...
            R_Topology_AllValues, ...
            R_Topology_NZValues, ...
            R_Weights_AllEdges, ...
            R_Weights_ConsensusEdges, ...
            Z_Topology_AllValues, ...
            Z_Topology_NZValues, ...
            Z_Weights_AllEdges, ...
            Z_Weights_ConsensusEdges, ...
            'VariableNames',{ ...
            'Outer', ...
            'MC', ...
            'R_Topology_AllValues', ...
            'R_Topology_NZValues', ...
            'R_Weights_AllEdges', ...
            'R_Weights_ConsensusEdges', ...
            'Z_Topology_AllValues', ...
            'Z_Topology_NZValues', ...
            'Z_Weights_AllEdges', ...
            'Z_Weights_ConsensusEdges'});

        RobustnessProfileSensitivity = ...
            [RobustnessProfileSensitivity;tmp];

    end
end



%% Statistical comparison of spectral-profile sensitivity

% Does the effect of changing Topology depend on the weight aggregation?
D_Topology = ...
    RobustnessProfileSensitivity.Z_Topology_AllValues - ...
    RobustnessProfileSensitivity.Z_Topology_NZValues;

[p_TopologyWeightDependence,~,stats_TopologyWeightDependence] = ...
    signrank(D_Topology,0);

% Does the effect of changing Weights depend on the topology?
D_Weights = ...
    RobustnessProfileSensitivity.Z_Weights_AllEdges - ...
    RobustnessProfileSensitivity.Z_Weights_ConsensusEdges;

[p_WeightTopologyDependence,~,stats_WeightTopologyDependence] = ...
    signrank(D_Weights,0);

% Average profile similarity following a Topology change
Z_Topology = mean([ ...
    RobustnessProfileSensitivity.Z_Topology_AllValues, ...
    RobustnessProfileSensitivity.Z_Topology_NZValues],2);

% Average profile similarity following a Weights change
Z_Weights = mean([ ...
    RobustnessProfileSensitivity.Z_Weights_AllEdges, ...
    RobustnessProfileSensitivity.Z_Weights_ConsensusEdges],2);

% Positive values indicate that changing Weights preserves the harmonic
% profile more strongly than changing Topology
D_Construction = Z_Weights - Z_Topology;

[p_Construction,~,stats_Construction] = ...
    signrank(D_Construction,0);



%% 2. Representativeness across construction schemes
%
% Examines the influence of connectome construction on Representativeness.
%
% Two construction factors are considered:
%
%   Topology:
%       AllEdges
%       ConsensusEdges
%
%   Weights:
%       AllValues
%       NZValues
%
% Three complementary analyses are performed:
%
%   1. Construction effects on Representativeness at the largest population
%      sample size (N = 700).
%
%   2. Construction effects on Representativeness gains across predefined
%      sample size contrasts.
%
%   3. Sensitivity of harmonic-wise Representativeness profiles to changes
%      in topology and weight aggregation at N = 700.


%% Endpoint construction effects

load(fullfile(ResultsDir,...
    'Representativeness_HarmonicSimilarityMatch_Results.mat'),...
    'RepresentativenessStats_AcrossHarmonics');

load(fullfile(ResultsDir,...
    'Representativeness_HarmonicSimilarityMatch_Results.mat'),'RepMeas');

load(fullfile(ResultsDir,...
    'Representativeness_HarmonicSimilarityMatch_Results.mat'),'RepMeasSubjects');

% Largest population sample size
Nmax_Rep = 700;

% Ensures construction factors are categorical
RepresentativenessStats_AcrossHarmonics.Topology = categorical( ...
    RepresentativenessStats_AcrossHarmonics.Topology);

RepresentativenessStats_AcrossHarmonics.Weights = categorical( ...
    RepresentativenessStats_AcrossHarmonics.Weights);

% Endpoint observations
T700 = RepresentativenessStats_AcrossHarmonics( ...
    RepresentativenessStats_AcrossHarmonics.N == Nmax_Rep,:);

% At N = 700, averages Monte Carlo occurrences while preserving
% Outer x Subject x construction
[G,TGroups] = findgroups( ...
    T700(:,{'Outer','Subject','Topology','Weights'}));

Median_HarmonicSimilarityMatch = splitapply( ...
    @mean, ...
    T700.Median_HarmonicSimilarityMatch, ...
    G);

% Number of Monte Carlo occurrences contributing to each value
NumMC = splitapply(@numel,T700.MC,G);

% Reduced endpoint table
RepresentativenessConstructionEndpointStats = [ ...
    TGroups, ...
    table(Median_HarmonicSimilarityMatch,NumMC)];

% Factorial mixed-effects model
RepresentativenessConstructionEndpointModel = fitlme( ...
    RepresentativenessConstructionEndpointStats, ...
    ['Median_HarmonicSimilarityMatch ~ Topology*Weights + ' ...
     '(1|Outer) + (1|Subject)']);

% Simple construction effects
R1 = GetLMEContrast( ...
    RepresentativenessConstructionEndpointModel,[0 1 0 0]);

R2 = GetLMEContrast( ...
    RepresentativenessConstructionEndpointModel,[0 1 0 1]);

R3 = GetLMEContrast( ...
    RepresentativenessConstructionEndpointModel,[0 0 1 0]);

R4 = GetLMEContrast( ...
    RepresentativenessConstructionEndpointModel,[0 0 1 1]);

% Combines endpoint construction effects
RepresentativenessConstructionEndpoint = [R1;R2;R3;R4];

RepresentativenessConstructionEndpoint.Effect = categorical({ ...
    'Topology_AllValues'; ...
    'Topology_NZValues'; ...
    'Weights_AllEdges'; ...
    'Weights_ConsensusEdges'});



%% Construction effects on sample size contrasts

% Representativeness outcome
Outcome = 'Median_HarmonicSimilarityMatch';

% Predefined sample size contrasts
NContrasts = [ ...
     10   50
     50  100
    100  200
    200  700];

% Will contain matched observations for all contrasts
RepresentativenessConstructionContrastStats = table();

% Constructs sample size gains
for nc = 1:size(NContrasts,1)

    % Current contrast
    N1 = NContrasts(nc,1);
    N2 = NContrasts(nc,2);

    % Observations at the two sample sizes
    T1raw = RepresentativenessStats_AcrossHarmonics( ...
        RepresentativenessStats_AcrossHarmonics.N == N1,:);

    T2raw = RepresentativenessStats_AcrossHarmonics( ...
        RepresentativenessStats_AcrossHarmonics.N == N2,:);


    % Average Monte Carlo occurrences while preserving
    % Outer x Subject x construction
    [G1,TGroups1] = findgroups( ...
        T1raw(:,{'Outer','Subject','Topology','Weights'}));

    Value_N1 = splitapply(@mean,T1raw.(Outcome),G1);

    T1 = [TGroups1,table(Value_N1)];


    [G2,TGroups2] = findgroups( ...
        T2raw(:,{'Outer','Subject','Topology','Weights'}));

    Value_N2 = splitapply(@mean,T2raw.(Outcome),G2);

    T2 = [TGroups2,table(Value_N2)];


    % Explicitly matches identical Outer x Subject x construction cases
    TD = innerjoin(T1,T2, ...
        'Keys',{'Outer','Subject','Topology','Weights'});

    % Ensures complete one-to-one matching
    assert(height(TD) == height(T1));
    assert(height(TD) == height(T2));

    % Larger harmonic similarity indicates greater Representativeness, so
    % positive values indicate improvement with increasing sample size
    TD.Gain = TD.Value_N2 - TD.Value_N1;

    % Stores contrast information
    TD.N1 = repmat(N1,height(TD),1);
    TD.N2 = repmat(N2,height(TD),1);
    TD.Contrast = repmat(nc,height(TD),1);

    % Appends observations
    RepresentativenessConstructionContrastStats = ...
        [RepresentativenessConstructionContrastStats;TD];
end

% Initializes model and result containers
RepresentativenessConstructionContrastModels = ...
    cell(size(NContrasts,1),1);

ContrastResults = ...
    cell(size(NContrasts,1),1);


% Fits one factorial model for each sample size contrast
for nc = 1:size(NContrasts,1)

    % Current contrast
    TD = RepresentativenessConstructionContrastStats( ...
        RepresentativenessConstructionContrastStats.Contrast == nc,:);

    % Factorial mixed-effects model
    mdl = fitlme(TD, ...
        ['Gain ~ Topology*Weights + ' ...
         '(1|Outer) + (1|Subject)']);

    RepresentativenessConstructionContrastModels{nc} = mdl;

    % Simple construction effects
    R1 = GetLMEContrast(mdl,[0 1 0 0]);
    R2 = GetLMEContrast(mdl,[0 1 0 1]);
    R3 = GetLMEContrast(mdl,[0 0 1 0]);
    R4 = GetLMEContrast(mdl,[0 0 1 1]);

    R = [R1;R2;R3;R4];

    % Labels construction effects
    R.Effect = categorical({ ...
        'Topology_AllValues'; ...
        'Topology_NZValues'; ...
        'Weights_AllEdges'; ...
        'Weights_ConsensusEdges'});

    % Stores contrast information
    R.N1 = repmat(NContrasts(nc,1),4,1);
    R.N2 = repmat(NContrasts(nc,2),4,1);
    R.Contrast = repmat(nc,4,1);

    ContrastResults{nc} = R;
end


% Combines results across sample size contrasts
RepresentativenessConstructionContrasts = ...
    vertcat(ContrastResults{:});

RepresentativenessConstructionContrasts = movevars( ...
    RepresentativenessConstructionContrasts, ...
    {'Contrast','N1','N2','Effect'}, ...
    'Before',1);



%% Spectral-profile sensitivity to construction scheme

% Finds N = 700 in the original Representativeness structure
SampleSizes_Rep = Cfg.Population.SampleSizes;
iN = find(SampleSizes_Rep == Nmax_Rep,1);

assert(~isempty(iN), ...
    'N = 700 was not found in the Representativeness sample sizes.');

% Will contain one observation per Outer x Subject occurrence
RepresentativenessProfileSensitivity = table();

for o = 1:numel(RepMeas)

    % Subject IDs represented in this outer split
    SubjectIDs = RepMeasSubjects{o}(:);

    % Number of left-out subjects
    nSubj = numel(SubjectIDs);


    % Harmonic-wise Representativeness profiles at N = 700, averaged
    % across Monte Carlo occurrences
    P_AE_AV = squeeze(mean( ...
        RepMeas{o}.AllEdges.AllValues(:,iN,:,:), ...
        1,'omitnan'));

    P_AE_NZ = squeeze(mean( ...
        RepMeas{o}.AllEdges.NZValues(:,iN,:,:), ...
        1,'omitnan'));

    P_CE_AV = squeeze(mean( ...
        RepMeas{o}.ConsensusEdges.AllValues(:,iN,:,:), ...
        1,'omitnan'));

    P_CE_NZ = squeeze(mean( ...
        RepMeas{o}.ConsensusEdges.NZValues(:,iN,:,:), ...
        1,'omitnan'));


    % Expected dimensions: subjects x harmonics
    assert(size(P_AE_AV,1) == nSubj);
    assert(size(P_AE_NZ,1) == nSubj);
    assert(size(P_CE_AV,1) == nSubj);
    assert(size(P_CE_NZ,1) == nSubj);

    % All construction schemes must contain the same number of harmonics
    assert(size(P_AE_AV,2) == size(P_AE_NZ,2));
    assert(size(P_AE_AV,2) == size(P_CE_AV,2));
    assert(size(P_AE_AV,2) == size(P_CE_NZ,2));


    % Construction sensitivity is assessed separately for each subject
    for is = 1:nSubj

        % Subject ID
        s = SubjectIDs(is);

        % Harmonic-wise profiles
        p_AE_AV = P_AE_AV(is,:)';
        p_AE_NZ = P_AE_NZ(is,:)';
        p_CE_AV = P_CE_AV(is,:)';
        p_CE_NZ = P_CE_NZ(is,:)';


        % Topology sensitivity:
        % changes Topology while keeping Weights fixed
        R_Topology_AllValues = corr( ...
            p_AE_AV,p_CE_AV, ...
            'Type','Pearson', ...
            'Rows','complete');

        R_Topology_NZValues = corr( ...
            p_AE_NZ,p_CE_NZ, ...
            'Type','Pearson', ...
            'Rows','complete');


        % Weight sensitivity:
        % changes Weights while keeping Topology fixed
        R_Weights_AllEdges = corr( ...
            p_AE_AV,p_AE_NZ, ...
            'Type','Pearson', ...
            'Rows','complete');

        R_Weights_ConsensusEdges = corr( ...
            p_CE_AV,p_CE_NZ, ...
            'Type','Pearson', ...
            'Rows','complete');


        % Avoids numerical problems if a correlation equals +/-1
        cliprho = @(r) min(max(r,-1+eps),1-eps);

        % Fisher Z transformation
        Z_Topology_AllValues = ...
            atanh(cliprho(R_Topology_AllValues));

        Z_Topology_NZValues = ...
            atanh(cliprho(R_Topology_NZValues));

        Z_Weights_AllEdges = ...
            atanh(cliprho(R_Weights_AllEdges));

        Z_Weights_ConsensusEdges = ...
            atanh(cliprho(R_Weights_ConsensusEdges));


        % Stores one observation per Outer x Subject occurrence
        tmp = table( ...
            o, ...
            s, ...
            R_Topology_AllValues, ...
            R_Topology_NZValues, ...
            R_Weights_AllEdges, ...
            R_Weights_ConsensusEdges, ...
            Z_Topology_AllValues, ...
            Z_Topology_NZValues, ...
            Z_Weights_AllEdges, ...
            Z_Weights_ConsensusEdges, ...
            'VariableNames',{ ...
            'Outer', ...
            'Subject', ...
            'R_Topology_AllValues', ...
            'R_Topology_NZValues', ...
            'R_Weights_AllEdges', ...
            'R_Weights_ConsensusEdges', ...
            'Z_Topology_AllValues', ...
            'Z_Topology_NZValues', ...
            'Z_Weights_AllEdges', ...
            'Z_Weights_ConsensusEdges'});

        RepresentativenessProfileSensitivity = ...
            [RepresentativenessProfileSensitivity;tmp];

    end
end



%% Subject-level aggregation of spectral-profile sensitivity

% Subjects can occur in multiple outer test sets. Repeated occurrences are
% therefore averaged for each subject in Fisher-Z space.
[G,Subject] = findgroups( ...
    RepresentativenessProfileSensitivity.Subject);

NumOccurrences = splitapply( ...
    @numel, ...
    RepresentativenessProfileSensitivity.Subject, ...
    G);

Z_Topology_AllValues = splitapply( ...
    @mean, ...
    RepresentativenessProfileSensitivity.Z_Topology_AllValues, ...
    G);

Z_Topology_NZValues = splitapply( ...
    @mean, ...
    RepresentativenessProfileSensitivity.Z_Topology_NZValues, ...
    G);

Z_Weights_AllEdges = splitapply( ...
    @mean, ...
    RepresentativenessProfileSensitivity.Z_Weights_AllEdges, ...
    G);

Z_Weights_ConsensusEdges = splitapply( ...
    @mean, ...
    RepresentativenessProfileSensitivity.Z_Weights_ConsensusEdges, ...
    G);


% Subject-level table
RepresentativenessProfileSensitivitySubject = table( ...
    Subject, ...
    NumOccurrences, ...
    Z_Topology_AllValues, ...
    Z_Topology_NZValues, ...
    Z_Weights_AllEdges, ...
    Z_Weights_ConsensusEdges);


% Corresponding correlations for descriptive interpretation
RepresentativenessProfileSensitivitySubject.R_Topology_AllValues = ...
    tanh(RepresentativenessProfileSensitivitySubject.Z_Topology_AllValues);

RepresentativenessProfileSensitivitySubject.R_Topology_NZValues = ...
    tanh(RepresentativenessProfileSensitivitySubject.Z_Topology_NZValues);

RepresentativenessProfileSensitivitySubject.R_Weights_AllEdges = ...
    tanh(RepresentativenessProfileSensitivitySubject.Z_Weights_AllEdges);

RepresentativenessProfileSensitivitySubject.R_Weights_ConsensusEdges = ...
    tanh(RepresentativenessProfileSensitivitySubject.Z_Weights_ConsensusEdges);



%% Statistical comparison of spectral-profile sensitivity

% Does the effect of changing Topology depend on the weight aggregation?
D_Topology = ...
    RepresentativenessProfileSensitivitySubject.Z_Topology_AllValues - ...
    RepresentativenessProfileSensitivitySubject.Z_Topology_NZValues;

[p_TopologyWeightDependence_Rep,~, ...
    stats_TopologyWeightDependence_Rep] = ...
    signrank(D_Topology,0);


% Does the effect of changing Weights depend on the topology?
D_Weights = ...
    RepresentativenessProfileSensitivitySubject.Z_Weights_AllEdges - ...
    RepresentativenessProfileSensitivitySubject.Z_Weights_ConsensusEdges;

[p_WeightTopologyDependence_Rep,~, ...
    stats_WeightTopologyDependence_Rep] = ...
    signrank(D_Weights,0);


% Average profile similarity following a Topology change
Z_Topology = mean([ ...
    RepresentativenessProfileSensitivitySubject.Z_Topology_AllValues, ...
    RepresentativenessProfileSensitivitySubject.Z_Topology_NZValues],2);

% Average profile similarity following a Weights change
Z_Weights = mean([ ...
    RepresentativenessProfileSensitivitySubject.Z_Weights_AllEdges, ...
    RepresentativenessProfileSensitivitySubject.Z_Weights_ConsensusEdges],2);

% Positive values indicate that changing Weights preserves the harmonic
% profile more strongly than changing Topology
D_Construction = Z_Weights - Z_Topology;

[p_Construction_Rep,~,stats_Construction_Rep] = ...
    signrank(D_Construction,0);



%% 3. Practical Utility across construction schemes
%
% Examines whether the frequency-dependent Practical Utility effect is
% preserved across connectome construction schemes.
%
% Four Practical Utility quantities are considered:
%
%   StructuralMismatch    Structural mismatch
%   Error_Functional      Functional error
%   Rel_DeltaC            Relative second-order moment error
%   Error_FC              FC error
%
% For each subject, metric and construction scheme, the difference between
% High- and Low-pass values is computed:
%
%   HighMinusLow = High - Low
%
% Repeated occurrences of the same subject across outer splits are first
% averaged. Construction effects are then assessed at the subject level.


%% Low-versus-High frequency effects

load(fullfile(ResultsDir,'PracticalUtility_Part1_Results.mat'),'PracticalUtilityStats');

% Ensures construction and filtering factors are categorical
PracticalUtilityStats.Topology = categorical( ...
    PracticalUtilityStats.Topology);

PracticalUtilityStats.Weights = categorical( ...
    PracticalUtilityStats.Weights);

PracticalUtilityStats.Filter = categorical( ...
    PracticalUtilityStats.Filter);

% Practical Utility outcomes
PUMetrics = { ...
    'StructuralMismatch', ...
    'Error_Functional', ...
    'Rel_DeltaC', ...
    'Error_FC'};

% Initializes output tables
PracticalUtilityConstructionStats = table();
PracticalUtilityConstructionSubject = table();
PracticalUtilityConstructionTests = table();


% Same analysis for each Practical Utility quantity
for m = 1:numel(PUMetrics)

    % Current outcome
    Outcome = PUMetrics{m};


    %% Low-versus-High frequency effect

    % Low-pass observations
    T_Low = PracticalUtilityStats( ...
        PracticalUtilityStats.Filter == "Low", ...
        {'Outer','Subject','Topology','Weights',Outcome});

    % High-pass observations
    T_High = PracticalUtilityStats( ...
        PracticalUtilityStats.Filter == "High", ...
        {'Outer','Subject','Topology','Weights',Outcome});

    % Renames outcome before matching
    T_Low.Properties.VariableNames{Outcome} = 'Low';
    T_High.Properties.VariableNames{Outcome} = 'High';

    % Explicitly matches Low and High observations within identical
    % Outer x Subject x construction cases
    TD = innerjoin(T_Low,T_High, ...
        'Keys',{'Outer','Subject','Topology','Weights'});

    % Ensures complete one-to-one matching
    assert(height(TD) == height(T_Low));
    assert(height(TD) == height(T_High));

    % High-minus-Low frequency effect
    TD.HighMinusLow = TD.High - TD.Low;


    %% Subject-level aggregation

    % Subjects can occur in multiple outer splits. Repeated occurrences are
    % averaged within each subject and construction scheme.
    [G,TGroups] = findgroups( ...
        TD(:,{'Subject','Topology','Weights'}));

    HighMinusLow = splitapply( ...
        @mean, ...
        TD.HighMinusLow, ...
        G);

    NumOccurrences = splitapply( ...
        @numel, ...
        TD.HighMinusLow, ...
        G);

    % Subject-level frequency effects
    TStats = [ ...
        TGroups, ...
        table(NumOccurrences,HighMinusLow)];

    % Metric identity
    TStats.Metric = repmat(string(Outcome),height(TStats),1);
    TStats = movevars(TStats,'Metric','Before',1);

    % Appends current metric
    PracticalUtilityConstructionStats = ...
        [PracticalUtilityConstructionStats;TStats];


    %% Construction effects

    % Converts the four construction schemes into paired subject-level
    % columns
    Subjects = unique(TStats.Subject);

    TSubject = table();

    for is = 1:numel(Subjects)

        % Current subject
        s = Subjects(is);

        % Corresponding construction schemes
        TS = TStats(TStats.Subject == s,:);

        % Extracts High-minus-Low effect for each construction scheme
        D_AE_AV = TS.HighMinusLow( ...
            TS.Topology == "AllEdges" & ...
            TS.Weights == "AllValues");

        D_AE_NZ = TS.HighMinusLow( ...
            TS.Topology == "AllEdges" & ...
            TS.Weights == "NZValues");

        D_CE_AV = TS.HighMinusLow( ...
            TS.Topology == "ConsensusEdges" & ...
            TS.Weights == "AllValues");

        D_CE_NZ = TS.HighMinusLow( ...
            TS.Topology == "ConsensusEdges" & ...
            TS.Weights == "NZValues");

        % Retains only subjects represented under all four construction
        % schemes
        if numel(D_AE_AV) ~= 1 || ...
           numel(D_AE_NZ) ~= 1 || ...
           numel(D_CE_AV) ~= 1 || ...
           numel(D_CE_NZ) ~= 1

            continue
        end

        % Stores paired construction values
        tmp = table( ...
            string(Outcome), ...
            s, ...
            D_AE_AV, ...
            D_AE_NZ, ...
            D_CE_AV, ...
            D_CE_NZ, ...
            'VariableNames',{ ...
            'Metric', ...
            'Subject', ...
            'D_AE_AV', ...
            'D_AE_NZ', ...
            'D_CE_AV', ...
            'D_CE_NZ'});

        TSubject = [TSubject;tmp];

    end


    % Effect of changing Topology while keeping Weights fixed
    TSubject.Topology_AllValues = ...
        TSubject.D_CE_AV - ...
        TSubject.D_AE_AV;

    TSubject.Topology_NZValues = ...
        TSubject.D_CE_NZ - ...
        TSubject.D_AE_NZ;

    % Effect of changing Weights while keeping Topology fixed
    TSubject.Weights_AllEdges = ...
        TSubject.D_AE_NZ - ...
        TSubject.D_AE_AV;

    TSubject.Weights_ConsensusEdges = ...
        TSubject.D_CE_NZ - ...
        TSubject.D_CE_AV;



    %% Relative construction effects

    % Expresses the magnitude of each construction-induced change relative
    % to the High-minus-Low effect for the corresponding reference
    % construction.

    % Topology changes
    TSubject.RelTopology_AllValues = ...
        abs(TSubject.Topology_AllValues) ./ ...
        abs(TSubject.D_AE_AV);

    TSubject.RelTopology_NZValues = ...
        abs(TSubject.Topology_NZValues) ./ ...
        abs(TSubject.D_AE_NZ);

    % Weight changes
    TSubject.RelWeights_AllEdges = ...
        abs(TSubject.Weights_AllEdges) ./ ...
        abs(TSubject.D_AE_AV);

    TSubject.RelWeights_ConsensusEdges = ...
        abs(TSubject.Weights_ConsensusEdges) ./ ...
        abs(TSubject.D_CE_AV);

    % Checks whether every subject had all four construction schemes
    assert(height(TSubject) == numel(Subjects), ...
        ['Some subjects did not have observations for all four ' ...
         'construction schemes for ',Outcome,'.']);

    % Appends current metric
    PracticalUtilityConstructionSubject = ...
        [PracticalUtilityConstructionSubject;TSubject];



    %% Statistical analyses

    % Tests whether changing Topology modifies the High-minus-Low effect
    [p_Topology_AllValues,~,stats_Topology_AllValues] = ...
        signrank(TSubject.Topology_AllValues,0);
    
    [p_Topology_NZValues,~,stats_Topology_NZValues] = ...
        signrank(TSubject.Topology_NZValues,0);
    
    % Tests whether changing Weights modifies the High-minus-Low effect
    [p_Weights_AllEdges,~,stats_Weights_AllEdges] = ...
        signrank(TSubject.Weights_AllEdges,0);
    
    [p_Weights_ConsensusEdges,~,stats_Weights_ConsensusEdges] = ...
        signrank(TSubject.Weights_ConsensusEdges,0);

    % Summarizes construction-effect tests
    TTests = table( ...
        repmat(string(Outcome),4,1), ...
        ["Topology_AllValues"; ...
         "Topology_NZValues"; ...
         "Weights_AllEdges"; ...
         "Weights_ConsensusEdges"], ...
        [median(TSubject.RelTopology_AllValues,'omitnan'); ...
         median(TSubject.RelTopology_NZValues,'omitnan'); ...
         median(TSubject.RelWeights_AllEdges,'omitnan'); ...
         median(TSubject.RelWeights_ConsensusEdges,'omitnan')], ...
        [stats_Topology_AllValues.zval; ...
         stats_Topology_NZValues.zval; ...
         stats_Weights_AllEdges.zval; ...
         stats_Weights_ConsensusEdges.zval], ...
        [p_Topology_AllValues; ...
         p_Topology_NZValues; ...
         p_Weights_AllEdges; ...
         p_Weights_ConsensusEdges], ...
        'VariableNames',{ ...
        'Metric', ...
        'Effect', ...
        'MedianRelativeChange', ...
        'Z', ...
        'pValue'});

    PracticalUtilityConstructionTests = ...
        [PracticalUtilityConstructionTests;TTests];
end



%% 4. Save results
% Saves all construction-scheme analysis outputs required for reporting
% and downstream figure generation

save(fullfile(ResultsDir,'ConstructionSchemeAnalysis_Results.mat'), ...
    'RobustnessConstructionEndpoint', ...
    'RobustnessConstructionContrastStats', ...
    'RobustnessConstructionContrasts', ...
    'RobustnessProfileSensitivity', ...
    'p_TopologyWeightDependence', ...
    'p_WeightTopologyDependence', ...
    'p_Construction', ...
    'RepresentativenessConstructionEndpointStats', ...
    'RepresentativenessConstructionEndpoint', ...
    'RepresentativenessConstructionContrastStats', ...
    'RepresentativenessConstructionContrasts', ...
    'RepresentativenessProfileSensitivity', ...
    'RepresentativenessProfileSensitivitySubject', ...
    'p_TopologyWeightDependence_Rep', ...
    'p_WeightTopologyDependence_Rep', ...
    'p_Construction_Rep', ...
    'PracticalUtilityConstructionStats', ...
    'PracticalUtilityConstructionSubject', ...
    'PracticalUtilityConstructionTests', ...
    '-v7.3');