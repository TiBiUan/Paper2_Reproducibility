function [Sparsity] = GSPTB_ComputeHarmonicDescriptors_Sparsity(u,T)
%GSPTB_COMPUTEHARMONICDESCRIPTORS_SPARSITY Compute the sparsity of one thresholded harmonic.
%
%   Sparsity = GSPTB_ComputeHarmonicDescriptors_Sparsity(u,T)
%
%   Description
%   -----------
%   Compute the sparsity of one thresholded harmonic.
%
%   Inputs
%   ------
%   u
%       R-element harmonic vector.
%   T
%       Non-negative absolute-amplitude threshold.
%
%   Outputs
%   -------
%   Sparsity
%       Fraction of harmonic entries set to zero after thresholding.
%
%   API status
%   ----------
%   Public.
%
%   See also
%   GSPTB_ComputeHarmonicDescriptors_Sparsity_sSC_Localization


    % We make the harmonic sparse according to threshold T
    u_sparse = u;
    u_sparse(abs(u)<T) = 0;

    % Sparsity
    Sparsity = sum(u_sparse==0)/length(u);
end
