function [DeltaSim,WithinSim,AcrossSim] = GSPTB_ComputeFamilyDeltaSim(CrossSubjSimMat,idx_fam)
%GSPTB_COMPUTEFAMILYDELTASIM Compare mean within-family and across-family harmonic similarity.
%
%   [DeltaSim,WithinSim,AcrossSim] = GSPTB_ComputeFamilyDeltaSim(CrossSubjSimMat,idx_fam)
%
%   Description
%   -----------
%   Compare mean within-family and across-family harmonic similarity.
%
%   Inputs
%   ------
%   CrossSubjSimMat
%       R-by-R harmonic similarity matrix.
%   idx_fam
%       R-element vector assigning each harmonic to a family indexed from 1 to F.
%
%   Outputs
%   -------
%   DeltaSim
%       Difference between the mean within-family and mean across-family
%       similarities.
%   WithinSim
%       F-element vector of within-family block means.
%   AcrossSim
%       F(F-1)/2-element vector of between-family block means.
%
%   API status
%   ----------
%   Public.
%
%   See also
%   GSPTB_ComputeFamilyDeltaNull, GSPTB_ExtractFamilies


    % Number of families
    F = max(idx_fam);

    % Initializes
    WithinSim = NaN(F,1);
    AcrossSim = NaN(F*(F-1)/2,1);

    % Within
    for f = 1:F
        tmp = CrossSubjSimMat(idx_fam==f,idx_fam==f);
        tmp = tmp(:);
        WithinSim(f) = mean(tmp);
    end
    
    idx_pairs = 1;
    
    for f1 = 1:F
        for f2 = 1:F
            if f2 > f1
                tmp = CrossSubjSimMat(idx_fam==f1,idx_fam==f2);
                tmp = tmp(:);
                AcrossSim(idx_pairs) = mean(tmp);
    
                idx_pairs = idx_pairs + 1;
            end
        end
    end
    
    DeltaSim = mean(WithinSim) - mean(AcrossSim);
end
