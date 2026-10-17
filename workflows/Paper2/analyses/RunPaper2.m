%% RunPaper2
%
% Reproduces the numerical analyses reported in Paper 2
%
% See README.md and DATA.md for installation and data requirements



%% 1. Configuration and setup

% Set the paths to the reproducibility repository and data archive
% EDIT THESE PATHS BEFORE RUNNING
RepositoryDir = '/path/to/Paper2_Reproducibility';
DataDir = '/path/to/Paper2_Data';

ResultsDir = fullfile(RepositoryDir,'results');

if ~isfolder(ResultsDir)
    mkdir(ResultsDir);
end

% Adds the required code
addpath(genpath(fullfile(RepositoryDir,'src')));
addpath(genpath(fullfile(RepositoryDir,'external')));
addpath(genpath(fullfile(RepositoryDir,'workflows','Paper2')));

% Gets configuration and local data paths
Cfg = GSPTB_Paper2Config();

% Loads the analysis-ready data
load(fullfile(DataDir,'Inputs','Data.mat'),'Data');

% Creates the population splits shared across analyses
load(fullfile(DataDir,'Inputs','PopSplits.mat'),'PopSplits');



%% 2. Robustness
Run_Robustness;



%% 3. Representativeness

% One call per measure to assess
RepCfg = struct;
RepCfg.MetricBaseName = 'HarmonicSimilarityMatch';
Run_Representativeness;

RepCfg.MetricBaseName = 'BandwidthIndex';
Run_Representativeness;

RepCfg.MetricBaseName = 'LeakageEntropy';
Run_Representativeness;

RepCfg.MetricBaseName = 'SubspaceSimilarityRequiredWindows';
Run_Representativeness;



%% 4. Practical Utility
Run_PracticalUtility_Part1;
Run_PracticalUtility_Part2;
Run_PracticalUtility_Part3;



%% 5. Cross-pillar analyses
Run_CrossPillarAnalysis;



%% 6. Connectome construction analyses
Run_ConstructionSchemeAnalysis;
