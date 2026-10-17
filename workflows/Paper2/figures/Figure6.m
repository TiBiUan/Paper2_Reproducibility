%% Figure 6 - Connectome construction
%
% Generates Figure 6 from previously computed analysis outputs.
%
% Figure 6A:
%   Harmonic-wise Robustness across connectome construction schemes.
%
% Figure 6B:
%   Harmonic-wise Representativeness across connectome construction schemes.
%
% Figure 6C:
%   Relative influence of connectome construction on the frequency-
%   dependent Practical Utility effect.
%
% This script performs visualization only. All analytical quantities and
% statistical comparisons are generated upstream.



%% Loading the required data

RepositoryDir = '/Users/bolton/Desktop/Paper2_Reproducibility';
DataDir = '/Users/bolton/Desktop/Paper2_Data';
ResultsDir = fullfile(RepositoryDir,'results');

addpath(genpath(fullfile(RepositoryDir,'src')));
addpath(genpath(fullfile(RepositoryDir,'external')));
addpath(genpath(fullfile(RepositoryDir,'workflows')));

Cfg = GSPTB_Paper2Config();

load(fullfile(ResultsDir,'Robustness_Results.mat'), ...
    'Robustness');

Rep = load(fullfile(ResultsDir,...
    'Representativeness_HarmonicSimilarityMatch_Results.mat'),'RepMeas');

load(fullfile(ResultsDir,'ConstructionSchemeAnalysis_Results.mat'), ...
    'PracticalUtilityConstructionSubject');


%% 2. Definitions

TopologyNames = ["AllEdges","ConsensusEdges"];
WeightsNames = ["AllValues","NZValues"];

SampleSizes_Rob = Cfg.Population.RobustnessSampleSizes;
SampleSizes_Rep = Cfg.Population.SampleSizes;

R = Cfg.Atlas.NumRegionsAnalysis;

% Practical Utility colors
cPU = [ ...
    102 194 165
    252 141  98
    141 160 203
    231 138 195] ./ 255;

% Heat-map color maps
CM_Reds = cbrewer2('seq','Reds',1000,'linear');
CM_Diff = cbrewer2('div','PiYG',1000,'linear');
CM_Greens = cbrewer2('seq','Greens',1000,'linear');


%% 3. Figure 6A - Robustness

% Harmonic x N x Outer x MC arrays
RobustnessHeatmaps = struct();

for t = 1:numel(TopologyNames)

    topology = TopologyNames(t);

    for w = 1:numel(WeightsNames)

        weights = WeightsNames(w);

        X = nan( ...
            R, ...
            numel(SampleSizes_Rob), ...
            numel(Robustness.Outer), ...
            numel(Robustness.Outer(1).MonteCarlo));

        for o = 1:numel(Robustness.Outer)

            for mc = 1:numel(Robustness.Outer(o).MonteCarlo)

                for iN = 1:numel(SampleSizes_Rob)

                    X(:,iN,o,mc) = ...
                        Robustness.Outer(o).MonteCarlo(mc). ...
                        SampleSize(iN). ...
                        (topology).(weights). ...
                        HarmonicSimilarity.Match(:);

                end
            end
        end

        % Displayed mean across realizations
        RobustnessHeatmaps.(topology).(weights) = ...
            mean(X,[3 4],'omitnan');

    end
end


% Construction differences
RobustnessHeatmaps.AllEdges.WeightDifference = ...
    RobustnessHeatmaps.AllEdges.AllValues - ...
    RobustnessHeatmaps.AllEdges.NZValues;

RobustnessHeatmaps.ConsensusEdges.WeightDifference = ...
    RobustnessHeatmaps.ConsensusEdges.AllValues - ...
    RobustnessHeatmaps.ConsensusEdges.NZValues;

RobustnessHeatmaps.TopologyDifference.AllValues = ...
    RobustnessHeatmaps.AllEdges.AllValues - ...
    RobustnessHeatmaps.ConsensusEdges.AllValues;

RobustnessHeatmaps.TopologyDifference.NZValues = ...
    RobustnessHeatmaps.AllEdges.NZValues - ...
    RobustnessHeatmaps.ConsensusEdges.NZValues;


%% Robustness heat maps

% Sample-size tick positions
TickN_Rob = [5 40 150 350];
[~,TickPos_Rob] = ismember(TickN_Rob,SampleSizes_Rob);

% Data to display
RobustnessPlotData = { ...
    RobustnessHeatmaps.AllEdges.AllValues, ...
    RobustnessHeatmaps.AllEdges.NZValues, ...
    RobustnessHeatmaps.AllEdges.WeightDifference, ...
    RobustnessHeatmaps.ConsensusEdges.AllValues, ...
    RobustnessHeatmaps.ConsensusEdges.NZValues, ...
    RobustnessHeatmaps.ConsensusEdges.WeightDifference, ...
    RobustnessHeatmaps.TopologyDifference.AllValues, ...
    RobustnessHeatmaps.TopologyDifference.NZValues};

RobustnessPlotTitles = { ...
    'All Values', ...
    'Non-zero Values', ...
    'All - Non-zero', ...
    '', ...
    '', ...
    '', ...
    '', ...
    ''};

% Absolute maps = red; difference maps = PiYG
RobustnessIsDifference = [false false true false false true true true];

for p = 1:numel(RobustnessPlotData)

    figure;

    imagesc( ...
        1:numel(SampleSizes_Rob), ...
        1:R, ...
        RobustnessPlotData{p});

    axis xy;

    % Sample-size axis
    xticks(TickPos_Rob);
    xticklabels(TickN_Rob);
    xlim([0.5 numel(SampleSizes_Rob)+0.5]);

    % Harmonic axis
    yticks([40 80 120 160 200]);
    ylim([0.5 R+0.5]);

    % Color scale
    if RobustnessIsDifference(p)

        colormap(gca,CM_Diff);
        clim([-0.5 0.5]);

    else

        colormap(gca,CM_Reds);
        clim([0 1]);

    end

    xlabel('Sample size');
    ylabel('Harmonics');

    if ~isempty(RobustnessPlotTitles{p})
        title(RobustnessPlotTitles{p});
    end

    box on;

end


%% 4. Figure 6B - Representativeness

RepMeas = Rep.RepMeas;

RepresentativenessHeatmaps = struct();

for t = 1:numel(TopologyNames)

    topology = TopologyNames(t);

    for w = 1:numel(WeightsNames)

        weights = WeightsNames(w);

        % Outer-specific harmonic x N means
        XOuter = nan( ...
            R, ...
            numel(SampleSizes_Rep), ...
            numel(RepMeas));

        for o = 1:numel(RepMeas)

            % MC x N x Subject x Harmonic
            X = RepMeas{o}.(topology).(weights);

            % Mean across MC and subjects
            X = mean(X,[1 3],'omitnan');

            % Harmonic x N
            X = squeeze(X)';

            XOuter(:,:,o) = X;

        end

        % Displayed mean across outer splits
        RepresentativenessHeatmaps.(topology).(weights) = ...
            mean(XOuter,3,'omitnan');

    end
end


% Construction differences
RepresentativenessHeatmaps.AllEdges.WeightDifference = ...
    RepresentativenessHeatmaps.AllEdges.AllValues - ...
    RepresentativenessHeatmaps.AllEdges.NZValues;

RepresentativenessHeatmaps.ConsensusEdges.WeightDifference = ...
    RepresentativenessHeatmaps.ConsensusEdges.AllValues - ...
    RepresentativenessHeatmaps.ConsensusEdges.NZValues;

RepresentativenessHeatmaps.TopologyDifference.AllValues = ...
    RepresentativenessHeatmaps.AllEdges.AllValues - ...
    RepresentativenessHeatmaps.ConsensusEdges.AllValues;

RepresentativenessHeatmaps.TopologyDifference.NZValues = ...
    RepresentativenessHeatmaps.AllEdges.NZValues - ...
    RepresentativenessHeatmaps.ConsensusEdges.NZValues;


%% Representativeness heat maps

% Sample-size tick positions
TickN_Rep = [5 40 150 350 600];
[~,TickPos_Rep] = ismember(TickN_Rep,SampleSizes_Rep);

% Data to display
RepresentativenessPlotData = { ...
    RepresentativenessHeatmaps.AllEdges.AllValues, ...
    RepresentativenessHeatmaps.AllEdges.NZValues, ...
    RepresentativenessHeatmaps.AllEdges.WeightDifference, ...
    RepresentativenessHeatmaps.ConsensusEdges.AllValues, ...
    RepresentativenessHeatmaps.ConsensusEdges.NZValues, ...
    RepresentativenessHeatmaps.ConsensusEdges.WeightDifference, ...
    RepresentativenessHeatmaps.TopologyDifference.AllValues, ...
    RepresentativenessHeatmaps.TopologyDifference.NZValues};

RepresentativenessPlotTitles = { ...
    'All Values', ...
    'Non-zero Values', ...
    'All - Non-zero', ...
    '', ...
    '', ...
    '', ...
    '', ...
    ''};

% Absolute maps = green; difference maps = PiYG
RepresentativenessIsDifference = ...
    [false false true false false true true true];

for p = 1:numel(RepresentativenessPlotData)

    figure;

    imagesc( ...
        1:numel(SampleSizes_Rep), ...
        1:R, ...
        RepresentativenessPlotData{p});

    axis xy;

    % Sample-size axis
    xticks(TickPos_Rep);
    xticklabels(TickN_Rep);
    xlim([0.5 numel(SampleSizes_Rep)+0.5]);

    % Harmonic axis
    yticks([40 80 120 160 200]);
    ylim([0.5 R+0.5]);

    % Color scale
    if RepresentativenessIsDifference(p)

        colormap(gca,CM_Diff);
        clim([-0.5 0.5]);

    else

        colormap(gca,CM_Greens);
        clim([0 1]);

    end

    xlabel('Sample size');
    ylabel('Harmonics');

    if ~isempty(RepresentativenessPlotTitles{p})
        title(RepresentativenessPlotTitles{p});
    end

    box on;

end


%% 5. Figure 6C - Practical Utility

T = PracticalUtilityConstructionSubject;

MetricNames = [ ...
    "StructuralMismatch", ...
    "Error_Functional", ...
    "Rel_DeltaC", ...
    "Error_FC"];

MetricLabels = { ...
    'Structural mismatch', ...
    'Functional error', ...
    '2nd-order moment error', ...
    'FC error'};

EffectNames = { ...
    'RelTopology_AllValues', ...
    'RelTopology_NZValues', ...
    'RelWeights_AllEdges', ...
    'RelWeights_ConsensusEdges'};

EffectLabels = { ...
    'Topology | All Values', ...
    'Topology | Non-zero Values', ...
    'Weights | All Edges', ...
    'Weights | Consensus Edges'};


figure;
hold on;

MetricPositions = [1 4 7 10];
WithinMetricSpacing = 0.45;

LegendHandles = gobjects(numel(EffectNames),1);


for m = 1:numel(MetricNames)

    Tm = T(string(T.Metric) == MetricNames(m),:);

    for e = 1:numel(EffectNames)

        % Relative construction effect, expressed as percentage
        y = 100 .* Tm.(EffectNames{e});

        y = y(isfinite(y));

        % Position within metric
        x0 = MetricPositions(m) + ...
            (e - 2.5) .* WithinMetricSpacing;

        % Density-dependent horizontal displacement
        [f,~] = ksdensity(y,y);

        if max(f) > 0
            f = f ./ max(f);
        end

        % Deterministic pseudo-jitter
        phi = (sqrt(5)-1)/2;
        u = mod((1:numel(y))' .* phi,1);

        x = x0 + ...
            (2*u-1) .* ...
            0.16 .* f(:);

        % Subject observations
        h = scatter( ...
            x, ...
            y, ...
            8, ...
            cPU(e,:), ...
            'filled', ...
            'MarkerFaceAlpha',0.20, ...
            'MarkerEdgeAlpha',0.20);

        % Distribution summary
        b = boxchart( ...
            repmat(x0,size(y)), ...
            y, ...
            'BoxFaceColor',cPU(e,:), ...
            'BoxEdgeColor',cPU(e,:), ...
            'WhiskerLineColor',cPU(e,:), ...
            'BoxWidth',0.30, ...
            'LineWidth',1.0, ...
            'MarkerStyle','none');

        if m == 1
            LegendHandles(e) = h;
        end

        uistack(b,'top');

    end
end


xlim([-0.5 11.5]);

xticks(MetricPositions);
xticklabels(MetricLabels);
xtickangle(0);

ylabel('Relative change in filtering effect [%]');

ylim([0,40]);

box off;

set(gca, ...
    'TickDir','out', ...
    'FontSize',12);

legend( ...
    LegendHandles, ...
    EffectLabels, ...
    'Box','on', ...
    'Location','northeast');