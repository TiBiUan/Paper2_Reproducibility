function [A_appr] = GSPTB_EstimateHarmonicRangeSC(U,Lambda,Deg,is_normalize)
%GSPTB_ESTIMATEHARMONICRANGESC Reconstruct structural connectivity from a selected set of harmonics.
%
%   A_appr = GSPTB_EstimateHarmonicRangeSC(U,Lambda,Deg,is_normalize)
%
%   Description
%   -----------
%   Reconstruct structural connectivity from a selected set of harmonics.
%
%   Inputs
%   ------
%   U
%       R-by-K matrix of selected harmonics.
%   Lambda
%       K-element vector of corresponding eigenvalues.
%   Deg
%       R-by-R degree matrix.
%   is_normalize
%       Logical flag controlling Frobenius normalization.
%
%   Outputs
%   -------
%   A_appr
%       R-by-R non-negative reconstructed structural connectivity matrix.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Negative reconstructed edge weights are set to zero.
%
%   See also
%   GSPTB_EstimateHarmonicSC, GSPTB_ComputeHarmonicRangeKernel


    % Rank-1 approximation of the normalized Laplacian using the harmonic
    % at hand
    L_appr = U*diag(Lambda)*U';

    DL = zeros(size(U,1),size(U,1));

    % Gets the diagonal of the approximated Laplacian
    for i = 1:size(U,1)
        DL(i,i) = L_appr(i,i);
    end

    % Computes the output
    A_appr = sqrt(Deg)*(DL-L_appr)*sqrt(Deg);

    % Sets negative-valued elements to zero
    A_appr(A_appr < 0) = 0;

    % Normalizes (Frobenius norm)
    if is_normalize
        A_appr = A_appr/norm(A_appr,'fro');
    end
end
