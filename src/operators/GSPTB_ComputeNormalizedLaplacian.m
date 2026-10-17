function [Ln] = GSPTB_ComputeNormalizedLaplacian(A)
%GSPTB_COMPUTENORMALIZEDLAPLACIAN Compute the symmetric normalized graph Laplacian.
%
%   Ln = GSPTB_ComputeNormalizedLaplacian(A)
%
%   Description
%   -----------
%   Compute the symmetric normalized graph Laplacian.
%
%   Inputs
%   ------
%   A
%       R-by-R weighted adjacency matrix.
%
%   Outputs
%   -------
%   Ln
%       R-by-R symmetric normalized Laplacian.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   The implementation assumes positive node degree for every node.
%
%   See also
%   GSPTB_ComputeSortedHarmonics


    k = sum(A,2);
    D = diag(k);

    % Standard Laplacian
    L = D - A;

    % Normalized Laplacian
    Ln = diag(1./sqrt(k)) * L * diag(1./sqrt(k));
end
