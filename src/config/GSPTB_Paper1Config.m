function [Cfg] = GSPTB_Paper1Config()

    % Acquisition
    Cfg.Data.TR = 0.72;
    Cfg.Data.ToLoad = 1:875;
    Cfg.Data.ScrubbingThreshold = 0.5;

    % Functional paradigms
    Cfg.Tasks.Names = { ...
        'EMOTION','GAMBLING','LANGUAGE','MOTOR', ...
        'RELATIONAL','SOCIAL','WM','REST1','REST2'};

    Cfg.Tasks.NumTimePoints = ...
        [170,247,310,278,226,268,399,1194,1194];

    Cfg.Tasks.MainFigureNames = { ...
        'EMOTION','GAMBLING','LANGUAGE','MOTOR', ...
        'RELATIONAL','SOCIAL','WM','REST1'};

    Cfg.Tasks.MainFigureIndices = ...
        find(ismember(Cfg.Tasks.Names,Cfg.Tasks.MainFigureNames));

    % Atlas
    Cfg.Atlas.IsCortexOnly = true;
    Cfg.Atlas.CorticalIndices = [1:108,136:243];
    Cfg.Atlas.Scale = 3;
    Cfg.Atlas.NumRegions = 216;

    % Structural connectomes
    Cfg.Structural.PropertyIndex = 3;

    % Harmonic analyses
    Cfg.Harmonics.NormalizeAdjacency = true;
    Cfg.Harmonics.ModularityGamma = 1;
    Cfg.Harmonics.SparsityThreshold = 0.001;
    Cfg.Harmonics.ReferenceSubject = 1;

    % Family organization
    Cfg.FamilyOrganization.ClusteringType = 'weighted';
    Cfg.FamilyOrganization.ClusteringDistance = 'cosine';
    Cfg.FamilyOrganization.CuttingCriterion = 'distance';
    Cfg.FamilyOrganization.CuttingDistance = 0.6;

    % Prediction
    Cfg.Prediction.Lambda = logspace(-4,4,400);
    Cfg.Prediction.NumOuterSplits = 100;
    Cfg.Prediction.NumInnerFolds = 3;
    Cfg.Prediction.LeftOutFraction = 0.2;

    % Task contrasts
    Cfg.TaskContrasts.TaskFactorName = 'Task';

    % Statistical inference
    Cfg.Statistics.Stability.NumNulls = 100;
    Cfg.Statistics.Stability.AlphaPercent = 5;
    Cfg.Statistics.Prediction.NumNulls = 1000;
    Cfg.Statistics.FamilyDelta.NumNulls = 100000;
    Cfg.Statistics.TaskContrasts.AlphaLevel = 0.05;
end