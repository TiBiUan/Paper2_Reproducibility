function [CrossSubjSimMat] = GSPTB_ComputeCrossSubjSimilarityMatrix(U)
%GSPTB_COMPUTECROSSSUBJSIMILARITYMATRIX Compute average all-to-all cross-subject harmonic similarity.
%
%   CrossSubjSimMat = GSPTB_ComputeCrossSubjSimilarityMatrix(U)
%
%   Description
%   -----------
%   Compute average all-to-all cross-subject harmonic similarity.
%
%   Inputs
%   ------
%   U
%       1-by-S cell array of aligned R-by-R harmonic bases.
%
%   Outputs
%   -------
%   CrossSubjSimMat
%       R-by-R matrix of mean absolute Pearson similarity across subject pairs for
%       every harmonic pair.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Absolute correlations reveal similarity between non-corresponding
%   harmonics beyond the matched diagonal.
%   For each harmonic pair, both subject-order orientations are pooled before
%   averaging.
%
%   Dependencies
%   ------------
%   jUpperTriMatToVec
%
%   See also
%   GSPTB_ComputeCrossSubjSimilarityCurve, GSPTB_ExtractFamilies


    % Number of regions
    R = size(U{1},1);

    % Number of subjects
    S = length(U);

    % Initialization
    CrossSubjSimMat = NaN(R,R);
    Harmo1 = NaN(R,S);
    Harmo2 = NaN(R,S);
    
    % Next, we expand the analyses to consider all pairs of harmonics
    for c1 = 1:R
    
        c1
    
        for s = 1:S
            Harmo1(:,s) = U{s}(:,c1);
        end
    
        for c2 = 1:R
    
            % We sample the two harmonics at play across subjects...
            for s = 1:S
                Harmo2(:,s) = U{s}(:,c2);
            end
        
            % ... and we compute the mean similarity across all pairs of
            % subjects, to fill our final matrix
            % Note: the absolute value is key to reveal the similarities beyond
            % the matching process!
            C = abs(corr(Harmo1,Harmo2,'Type','Pearson'));
            CrossSubjSimMat(c1,c2) = mean([jUpperTriMatToVec(C);jUpperTriMatToVec(C')]);
        end
    end
end
