function [PredSplits] = GSPTB_CreateBehavioralPredictionSplits(Opts)

% GSPTB_CreatePredictionSplits
%
% Creates mutually exclusive, family-aware folds for cross-fitted
% generation of population-derived features.
%
% Every subject:
%   - belongs to exactly one test fold;
%   - belongs to the training set of all other folds.
%
% Complete family components are always assigned to the same fold.
%
% Required inputs:
%   Opts.FamilyKernel : S x S family-relationship matrix
%   Opts.NumFolds     : number of cross-fitting folds
%
% Output:
%   PredSplits.Fold(k).Train
%   PredSplits.Fold(k).Test
%   PredSplits.Fold(k).NumTrain
%   PredSplits.Fold(k).NumTest
%   PredSplits.FoldID
%   PredSplits.FoldSize
%   PredSplits.NumFolds


%% Basic inputs

% Number of subjects
S = size(Opts.FamilyKernel,1);

% Number of desired folds
K = Opts.NumFolds;


%% Input checks

if size(Opts.FamilyKernel,2) ~= S
    error('Opts.FamilyKernel must be a square matrix.');
end

if K < 2 || mod(K,1) ~= 0
    error('Opts.NumFolds must be an integer >= 2.');
end

if K > S
    error('Opts.NumFolds cannot exceed the number of subjects.');
end


%% Identify complete family components

% Family relationships are represented as a graph
G = graph(Opts.FamilyKernel);

% Subjects within the same connected component belong to the same family
family = conncomp(G);

% Unique family identifiers
fam_ids = unique(family);

% Number of independent family components
nFamilies = numel(fam_ids);

if K > nFamilies
    error(['Opts.NumFolds cannot exceed the number of independent ' ...
        'family components.']);
end


%% Compute family sizes

FamSize = zeros(nFamilies,1);

for f = 1:nFamilies
    FamSize(f) = sum(family == fam_ids(f));
end


%% Order families for balanced allocation

% Randomize first so that families of equal size do not always receive
% the same deterministic ordering
perm = randperm(nFamilies);

fam_ids = fam_ids(perm);
FamSize = FamSize(perm);

% Assign large families first; this improves fold-size balancing
[~,order] = sort(FamSize,'descend');

fam_ids = fam_ids(order);
FamSize = FamSize(order);


%% Assign each complete family to exactly one fold

% Number of subjects currently assigned to every fold
FoldSize = zeros(1,K);

% Fold membership of every subject
FoldID = zeros(S,1);

for f = 1:nFamilies

    % All subjects belonging to the current family
    idx_family = find(family == fam_ids(f));

    % Select one of the currently smallest folds.
    %
    % If several folds have the same current size, select one randomly
    % to avoid systematic preference for low-numbered folds.
    min_size = min(FoldSize);
    candidate_folds = find(FoldSize == min_size);

    k = candidate_folds(randi(numel(candidate_folds)));

    % Assign the complete family to this fold
    FoldID(idx_family) = k;

    % Update fold size
    FoldSize(k) = FoldSize(k) + numel(idx_family);
end


%% Create train/test indices

for k = 1:K

    % Every subject assigned to this fold is test
    PredSplits.Fold(k).Test = find(FoldID == k);

    % Every other subject is training
    PredSplits.Fold(k).Train = find(FoldID ~= k);

    % Exact numbers
    PredSplits.Fold(k).NumTest = ...
        numel(PredSplits.Fold(k).Test);

    PredSplits.Fold(k).NumTrain = ...
        numel(PredSplits.Fold(k).Train);
end


%% Save metadata

PredSplits.NumFolds = K;
PredSplits.FoldID = FoldID;
PredSplits.FoldSize = FoldSize;

PredSplits.NumSubjects = S;
PredSplits.NumFamilies = nFamilies;


%% Sanity checks

% Collect all test indices
AllTest = [];

for k = 1:K
    AllTest = [AllTest; PredSplits.Fold(k).Test(:)];
end

% Every subject must appear in a test fold
assert(numel(AllTest) == S,...
    'Not every subject was assigned to exactly one test fold.');

% No subject may appear more than once
assert(numel(unique(AllTest)) == S,...
    'At least one subject appears in more than one test fold.');

% Test folds must collectively contain subjects 1:S
assert(isequal(sort(AllTest(:)),(1:S)'),...
    'Test folds do not reproduce the complete subject population.');

% No family may be split across folds
for f = 1:nFamilies

    idx_family = find(family == fam_ids(f));

    assert(numel(unique(FoldID(idx_family))) == 1,...
        'At least one family was split across folds.');
end


%% Display fold sizes

disp('Prediction cross-fitting fold sizes:');

for k = 1:K
    disp(['Fold ',num2str(k),': ',...
        num2str(PredSplits.Fold(k).NumTrain),' train / ',...
        num2str(PredSplits.Fold(k).NumTest),' test']);
end

end