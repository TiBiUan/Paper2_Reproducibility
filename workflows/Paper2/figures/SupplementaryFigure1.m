%% Loading the required data

RepositoryDir = '/Users/bolton/Desktop/Paper2_Reproducibility';
DataDir = '/Users/bolton/Desktop/Paper2_Data';
ResultsDir = fullfile(RepositoryDir,'results');

addpath(genpath(fullfile(RepositoryDir,'src')));
addpath(genpath(fullfile(RepositoryDir,'external')));

load(fullfile(ResultsDir,'Robustness_Results.mat'));



%% Supplementary Figure 1 - Robustness convergence

% Colors
CM_Reds = cbrewer2('seq','Reds',9,'linear');
CM_Oranges = cbrewer2('seq','Oranges',9,'linear');

CM_Gains = [ ...
    CM_Reds(6,:); ...      % Harmonic similarity
    CM_Oranges(5,:); ...   % Projection-error required windows
    CM_Reds(9,:)];         % Laplacian reconstruction error


%% Supplementary Figure 1A - Marginal gain rates

xdata = Nmid(:)';

for m = 1:3

    figure;
    hold on;

    % Individual {Outer,MC} trajectories
    for i = 1:size(GainRate{m},1)

        plot( ...
            xdata, ...
            GainRate{m}(i,:), ...
            'LineWidth',0.25, ...
            'Color',[0.6 0.6 0.6]);

    end

    % Ensures plotting quantities are row vectors
    y025 = GainP025{m}(:)';
    y975 = GainP975{m}(:)';
    yMed = GainMed{m}(:)';

    % 2.5th-97.5th percentile interval
    fill( ...
        [xdata fliplr(xdata)], ...
        [y025 fliplr(y975)], ...
        CM_Gains(m,:), ...
        'EdgeColor','none', ...
        'FaceAlpha',0.5);

    % Median across {Outer,MC} trajectories
    plot( ...
        xdata, ...
        yMed, ...
        'LineWidth',2, ...
        'Color',CM_Gains(m,:), ...
        'Marker','o', ...
        'MarkerFaceColor',CM_Gains(m,:), ...
        'MarkerSize',6, ...
        'MarkerEdgeColor','k');

    xlim([xdata(1) xdata(end)]);
    box off;
end


%% Supplementary Figure 1B - Power-law parameters

% -------------------------------------------------------------------------
% Convergence exponent p
% -------------------------------------------------------------------------

AlphaVal = 0.3;

figure;
hold on;

histogram( ...
    PowerLawParams.p(1,:), ...
    0.3:0.005:0.8, ...
    'FaceColor',CM_Gains(1,:), ...
    'FaceAlpha',AlphaVal);

histogram( ...
    PowerLawParams.p(2,:), ...
    0.3:0.005:0.8, ...
    'FaceColor',CM_Gains(2,:), ...
    'FaceAlpha',AlphaVal);

histogram( ...
    PowerLawParams.p(3,:), ...
    0.3:0.005:0.8, ...
    'FaceColor',CM_Gains(3,:), ...
    'FaceAlpha',AlphaVal);

xlabel('Parameter p');
ylabel('Counts');

box off;


% -------------------------------------------------------------------------
% Scaling parameter a
% -------------------------------------------------------------------------

figure;

for m = 1:3

    subplot(1,3,m);

    histogram( ...
        PowerLawParams.a(m,:), ...
        30, ...
        'FaceColor',CM_Gains(m,:));

    xlabel('Parameter a');
    ylabel('Counts');

    box off;
end