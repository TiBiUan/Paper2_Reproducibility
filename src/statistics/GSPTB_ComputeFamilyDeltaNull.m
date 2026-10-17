function [DeltaSimNull] = GSPTB_ComputeFamilyDeltaNull(CrossSubjSimMat,FamSize,n_nulls_famdelta)
%GSPTB_COMPUTEFAMILYDELTANULL Generate a null distribution for the within-minus-across-family similarity statistic.
%
%   DeltaSimNull = GSPTB_ComputeFamilyDeltaNull(CrossSubjSimMat,FamSize,n_nulls_famdelta)
%
%   Description
%   -----------
%   Generate a null distribution for the within-minus-across-family similarity
%   statistic.
%
%   Inputs
%   ------
%   CrossSubjSimMat
%       R-by-R harmonic similarity matrix.
%   FamSize
%       F-element vector of family sizes.
%   n_nulls_famdelta
%       Number of null realizations.
%
%   Outputs
%   -------
%   DeltaSimNull
%       n_nulls_famdelta-by-1 vector of null within-minus-across-family
%       differences.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Each null realization preserves family sizes and spectral contiguity while
%   randomizing the ordering of blocks.
%
%   Dependencies
%   ------------
%   GSPTB_ComputeFamilyDeltaSim
%
%   See also
%   GSPTB_ComputeFamilyDeltaSim


    % Number of families
    F = length(FamSize);

    DeltaSimNull = NaN(n_nulls_famdelta,1);

    % Null realizations
    for n = 1:n_nulls_famdelta
        
        % Shuffles the ordering of the blocks
        current_sizes = FamSize(randperm(F));
        
        id = 1;
    
        for f = 1:F
            current_idx(id:id+current_sizes(f)-1) = f;
            id = id + current_sizes(f);
        end
    
        DeltaSimNull(n) = GSPTB_ComputeFamilyDeltaSim(CrossSubjSimMat,current_idx);
    end
end
