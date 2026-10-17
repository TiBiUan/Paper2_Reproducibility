function [Ker_SC] = GSPTB_ComputeHarmonicRangeKernel(U,Lambda,DegMat,...
    is_normalize,Harm_Range)
%GSPTB_COMPUTEHARMONICRANGEKERNEL Construct a subject kernel from connectivity reconstructed over a harmonic range.
%
%   Ker_SC = GSPTB_ComputeHarmonicRangeKernel(U,Lambda,DegMat,is_normalize,Harm_Range)
%
%   Description
%   -----------
%   Construct a subject kernel from connectivity reconstructed over a harmonic
%   range.
%
%   Inputs
%   ------
%   U
%       1-by-S cell array of R-by-R harmonic bases.
%   Lambda
%       1-by-S cell array of R-element eigenvalue vectors.
%   DegMat
%       1-by-S cell array of R-by-R degree matrices.
%   is_normalize
%       Logical flag controlling Frobenius normalization of reconstructed
%       connectivity.
%   Harm_Range
%       Vector of harmonic indices retained in the reconstruction.
%
%   Outputs
%   -------
%   Ker_SC
%       S-by-S Pearson-correlation kernel between vectorized reconstructed
%       connectivity matrices.
%
%   API status
%   ----------
%   Public.
%
%   Dependencies
%   ------------
%   GSPTB_EstimateHarmonicRangeSC
%   jUpperTriMatToVec
%
%   See also
%   GSPTB_KernelRidgeRegression


    % Number of harmonics
    R = size(U{1},1);

    % Number of subjects
    S = length(U);

    % Initializes
    SC_harm = NaN(R*(R-1)/2,S);

    % Same for each subject...
    for s = 1:S

        % This computes the approximated SC matrix when using the
        % specified range of harmonics
        A_appr = GSPTB_EstimateHarmonicRangeSC(U{s}(:,Harm_Range),...
            Lambda{s}(Harm_Range),DegMat{s},is_normalize);
        
        % Gets the vectorized representations
        SC_harm(:,s) = jUpperTriMatToVec(A_appr);
    end

    % Gets the associated kernel
    Ker_SC = corr(SC_harm);
end
