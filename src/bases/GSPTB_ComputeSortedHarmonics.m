function [U_out,Lambda_out] = GSPTB_ComputeSortedHarmonics(Ln)
%GSPTB_COMPUTESORTEDHARMONICS Compute and sort the eigenvectors and eigenvalues of a graph operator.
%
%   [U_out,Lambda_out] = GSPTB_ComputeSortedHarmonics(Ln)
%
%   Description
%   -----------
%   Compute and sort the eigenvectors and eigenvalues of a graph operator.
%
%   Inputs
%   ------
%   Ln
%       R-by-R symmetric graph operator, typically a normalized Laplacian.
%
%   Outputs
%   -------
%   U_out
%       R-by-R matrix whose columns are eigenvectors ordered by ascending
%       eigenvalue.
%   Lambda_out
%       R-element vector of sorted eigenvalues.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   The routine uses a full eigendecomposition and assumes a real symmetric
%   input.
%
%   See also
%   GSPTB_ComputeNormalizedLaplacian


    [U,Lambda] = eig(Ln);
    [Lambda,idx] = sort(diag(Lambda),'ascend');
    U_out = U(:,idx);
    Lambda_out = Lambda;
end
