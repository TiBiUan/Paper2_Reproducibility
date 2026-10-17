function [A_appr] = GSPTB_EstimateHarmonicSC(u,lam,Deg,is_normalize)
%GSPTB_ESTIMATEHARMONICSC Reconstruct structural connectivity from a single harmonic.
%
%   A_appr = GSPTB_EstimateHarmonicSC(u,lam,Deg,is_normalize)
%
%   Description
%   -----------
%   Reconstruct structural connectivity from a single harmonic.
%
%   Inputs
%   ------
%   u
%       R-element harmonic vector.
%   lam
%       Associated scalar eigenvalue.
%   Deg
%       R-by-R degree matrix.
%   is_normalize
%       Logical flag controlling Frobenius normalization.
%
%   Outputs
%   -------
%   A_appr
%       R-by-R non-negative harmonic-wise reconstructed structural connectivity
%       matrix.
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
%   GSPTB_EstimateHarmonicRangeSC


    % Rank-1 approximation of the normalized Laplacian using the harmonic
    % at hand
    L_appr = lam*(u*u');

    DL = zeros(length(u),length(u));

    % Gets the diagonal of the approximated Laplacian
    for i = 1:length(u)
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
