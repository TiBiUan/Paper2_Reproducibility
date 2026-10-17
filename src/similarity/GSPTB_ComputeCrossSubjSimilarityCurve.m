function [CrossSubjSimCurve] = GSPTB_ComputeCrossSubjSimilarityCurve(U)
%GSPTB_COMPUTECROSSSUBJSIMILARITYCURVE Compute cross-subject similarity between matched harmonics.
%
%   CrossSubjSimCurve = GSPTB_ComputeCrossSubjSimilarityCurve(U)
%
%   Description
%   -----------
%   Compute cross-subject similarity between matched harmonics.
%
%   Inputs
%   ------
%   U
%       1-by-S cell array of aligned R-by-R harmonic bases.
%
%   Outputs
%   -------
%   CrossSubjSimCurve
%       S(S-1)/2-by-R matrix of pairwise Pearson correlations for each matched
%       harmonic.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Columns index harmonic position; rows index unique subject pairs.
%
%   Dependencies
%   ------------
%   jUpperTriMatToVec
%
%   See also
%   GSPTB_ComputeCrossSubjSimilarityMatrix


    % Number of regions
    R = size(U{1},1);

    % Number of subjects
    S = length(U);

    % Initialization
    CrossSubjSimCurve = NaN(S*(S-1)/2,R);
    tmp_CH = NaN(R,S);
    
    % First, we consider only matched harmonics
    for ch = 1:R
    
        % We fill a matrix (R x S) with the data: the same harmonic is taken 
        % across all subjects...
        for s = 1:S
            tmp_CH(:,s) = U{s}(:,ch);
        end
    
        % ... and then, we assess similarity by looking at all pairs of
        % subjects at once in a vector form
        CrossSubjSimCurve(:,ch) = jUpperTriMatToVec(corr(tmp_CH,'Type','Pearson'));
    end
end
