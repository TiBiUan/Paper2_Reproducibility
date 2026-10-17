%% Loading the required data

RepositoryDir = '/Users/bolton/Desktop/Paper2_Reproducibility';
DataDir = '/Users/bolton/Desktop/Paper2_Data';
ResultsDir = fullfile(RepositoryDir,'results');

addpath(genpath(fullfile(RepositoryDir,'src')));
addpath(genpath(fullfile(RepositoryDir,'external')));



%% Supplementary Figure 3 - Additional Representativeness measures

% Color map
CM_YlGn = cbrewer2('seq','YlGn',1000,'linear');
CM_Greens = cbrewer2('seq','Greens',1000,'linear');


%% Leakage entropy

load(fullfile(ResultsDir,'Representativeness_LeakageEntropy_Results.mat'), ...
    'RepresentativenessHarmonicDescriptives');

T = RepresentativenessHarmonicDescriptives( ...
    RepresentativenessHarmonicDescriptives.Topology == "AllEdges" & ...
    RepresentativenessHarmonicDescriptives.Weights == "AllValues",:);

NSamples = unique(T.N);
Harmonics = unique(T.Harmonic);

MeanMeasure = NaN(numel(Harmonics),numel(NSamples));
Stability = NaN(numel(Harmonics),numel(NSamples));

for n = 1:numel(NSamples)

    TN = T(T.N == NSamples(n),:);
    TN = sortrows(TN,'Harmonic');

    MeanMeasure(:,n) = TN.E_Realization_Mean_Subject;
    Stability(:,n) = TN.StabilityRatio;

end

% Mean measure value
figure;
imagesc(MeanMeasure);

colormap(CM_Greens);
clim([0 3.5]);

box off;


% Stability of inter-individual differences
figure;
imagesc(Stability);

colormap(CM_YlGn);
clim([0 1]);

box off;


%% Bandwidth index

load(fullfile(ResultsDir,'Representativeness_BandwidthIndex_Results.mat'), ...
    'RepresentativenessHarmonicDescriptives');

T = RepresentativenessHarmonicDescriptives( ...
    RepresentativenessHarmonicDescriptives.Topology == "AllEdges" & ...
    RepresentativenessHarmonicDescriptives.Weights == "AllValues",:);

NSamples = unique(T.N);
Harmonics = unique(T.Harmonic);

MeanMeasure = NaN(numel(Harmonics),numel(NSamples));
Stability = NaN(numel(Harmonics),numel(NSamples));

for n = 1:numel(NSamples)

    TN = T(T.N == NSamples(n),:);
    TN = sortrows(TN,'Harmonic');

    MeanMeasure(:,n) = TN.E_Realization_Mean_Subject;
    Stability(:,n) = TN.StabilityRatio;

end

% Mean measure value
figure;
imagesc(MeanMeasure);

colormap(CM_Greens);
clim([0 22]);

box off;


% Stability of inter-individual differences
figure;
imagesc(Stability);

colormap(CM_YlGn);
clim([0 1]);

box off;


%% Required subspace window

load(fullfile(ResultsDir,'Representativeness_SubspaceSimilarityRequiredWindows_Results.mat'), ...
    'RepresentativenessHarmonicDescriptives');

T = RepresentativenessHarmonicDescriptives( ...
    RepresentativenessHarmonicDescriptives.Topology == "AllEdges" & ...
    RepresentativenessHarmonicDescriptives.Weights == "AllValues",:);

NSamples = unique(T.N);
Harmonics = unique(T.Harmonic);

MeanMeasure = NaN(numel(Harmonics),numel(NSamples));
Stability = NaN(numel(Harmonics),numel(NSamples));

for n = 1:numel(NSamples)

    TN = T(T.N == NSamples(n),:);
    TN = sortrows(TN,'Harmonic');

    MeanMeasure(:,n) = TN.E_Realization_Mean_Subject;
    Stability(:,n) = TN.StabilityRatio;

end

% Mean measure value
figure;
imagesc(MeanMeasure);

colormap(CM_Greens);
clim([0 100]);

box off;


% Stability of inter-individual differences
figure;
imagesc(Stability);

colormap(CM_YlGn);
clim([0 1]);

box off;