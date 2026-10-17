function [Cfg] = GSPTB_Paper2Config()

    % Acquisition
    Cfg.Data.TR = 0.72;
    Cfg.Data.ToLoad = 1:875;
    Cfg.Data.ScrubbingThreshold = 0.5;

    % Resting-state runs
    Cfg.Rest.Names = {'REST1','REST2'};
    Cfg.Rest.NumTimePoints = [1194,1194];

    % Atlas
    Cfg.Atlas.IsCortexOnly = true;
    Cfg.Atlas.CorticalIndices = [1:108,136:243];
    Cfg.Atlas.Scale = 3;
    Cfg.Atlas.NumRegions = 274;
    Cfg.Atlas.NumRegionsAnalysis = 216;

    % Structural connectomes
    Cfg.Structural.PropertyIndex = 3;

    % Harmonic analyses
    Cfg.Harmonics.ReferenceSubject = 1;

    Cfg.Population.LeftOutFraction = 0.2;
    Cfg.Population.TrainSize = floor((1-Cfg.Population.LeftOutFraction)*numel(Cfg.Data.ToLoad));
    Cfg.Population.NumOuterSplits = 10;
    Cfg.Population.NumMonteCarlo = 20;
    Cfg.Population.SampleSizes = [5,10,20,30,40,50,75,100,125,150,175,200,250,300,350,400,450,500,550,600,650,700];
    Cfg.Population.RobustnessSampleSizes = ...
        Cfg.Population.SampleSizes(2*Cfg.Population.SampleSizes <= Cfg.Population.TrainSize);
    Cfg.Population.NumBins = 20;
    Cfg.Population.TopologySchemes = {'all','consensus'};
    Cfg.Population.WeightSchemes = {'all','nonzero'};
    
    Cfg.Population.ProjectionWindows = 0:Cfg.Atlas.NumRegionsAnalysis-1;
    Cfg.Population.ProjectionTolerance = 0.1;

    Cfg.Population.SubspaceWindows = [0:10,12:2:30,35:5:60,70:10:150];
    Cfg.Population.SubspaceSimilarityThreshold = 0.9;

    % Contrasts considered for sample sizes, based on recurring literature 
    % values
    Cfg.Population.RobustnessContrasts = [10,50; 50,100;100,200;200,350];
    Cfg.Population.RepresentativenessContrasts = [10,50; 50,100;100,200;200,350];
    Cfg.Population.PracticalUtilityContrasts = [10,50; 50,100;100,200;200,350];

    % Tangibility thresholds for Robustness analyses
    % HWSM and projection-window thresholds are expressed in native metric
    % units; the Laplacian reconstruction error threshold is a relative
    % percentage reduction
    Cfg.Robustness.Tangibility.Median_HarmSim_Match = 0.01;
    Cfg.Robustness.Tangibility.P95_ProjErr_ReqWin = 1;
    Cfg.Robustness.Tangibility.LaplacianReconstructionError = 5;

    % Tangibility threshold for frequency-dependent differences in the
    % primary normalized Practical Utility measures
    Cfg.PracticalUtility.TangibilityThreshold = 0.05;

    % Prediction
    Cfg.Prediction.Lambda = logspace(-4,4,400);
    Cfg.Prediction.NumInnerFolds = 3;
end