%% Loading the required data

RepositoryDir = '/Users/bolton/Desktop/Paper2_Reproducibility';
DataDir = '/Users/bolton/Desktop/Paper2_Data';
ResultsDir = fullfile(RepositoryDir,'results');

addpath(genpath(fullfile(RepositoryDir,'src')));
addpath(genpath(fullfile(RepositoryDir,'external')));

load(fullfile(ResultsDir,'Representativeness_HarmonicSimilarityMatch_Results.mat'));



%% Supplementary Figure 2 - Representativeness variability

% Reference connectome construction
T = RepresentativenessHarmonicDescriptives( ...
    RepresentativenessHarmonicDescriptives.Topology == "AllEdges" & ...
    RepresentativenessHarmonicDescriptives.Weights == "AllValues",:);

% Sample sizes and harmonics
NSamples = unique(T.N);
Harmonics = unique(T.Harmonic);

% Color maps
CM_YlGn = cbrewer2('seq','YlGn',1000,'linear');
CM_SampleSizes = cbrewer2( ...
    'div','PuOr',numel(NSamples),'linear');


%% Supplementary Figure 2A - Persistent between-subject variability

% Converts the plotting quantity from long to Harmonic x N form
Persistent = NaN(numel(Harmonics),numel(NSamples));

for n = 1:numel(NSamples)

    TN = T(T.N == NSamples(n),:);
    TN = sortrows(TN,'Harmonic');

    Persistent(:,n) = TN.SD_Subject_E_Realization;

end

figure;

imagesc(Persistent);

colormap(CM_YlGn);
clim([0 0.16]);

box off;


%% Supplementary Figure 2B - Typical within-realization variability

% Converts the plotting quantity from long to Harmonic x N form
Typical = NaN(numel(Harmonics),numel(NSamples));

for n = 1:numel(NSamples)

    TN = T(T.N == NSamples(n),:);
    TN = sortrows(TN,'Harmonic');

    Typical(:,n) = TN.E_Realization_SD_Subject;

end

figure;

imagesc(Typical);

colormap(CM_YlGn);
clim([0 0.16]);

box off;


%% Supplementary Figure 2C - Relationship between variability quantities

figure;
hold on;

% Identity line
plot([0 0.16],[0 0.16],'k--');

% One point per harmonic, colored according to sample size
for n = 1:numel(NSamples)

    TN = T(T.N == NSamples(n),:);

    scatter( ...
        TN.E_Realization_SD_Subject, ...
        TN.SD_Subject_E_Realization, ...
        5, ...
        CM_SampleSizes(n,:), ...
        'filled');

end

xlim([0 0.16]);
ylim([0 0.16]);

box off;