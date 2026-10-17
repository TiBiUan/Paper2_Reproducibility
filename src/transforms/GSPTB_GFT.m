function [Xhat] = GSPTB_GFT(U,X)
%GSPTB_GFT Compute the graph Fourier transform of graph signals.
%
%   Xhat = GSPTB_GFT(U,X)
%
%   Description
%   -----------
%   Compute the graph Fourier transform of graph signals.
%
%   Inputs
%   ------
%   U
%       R-by-R orthonormal graph spectral basis.
%   X
%       R-by-T graph-signal matrix.
%
%   Outputs
%   -------
%   Xhat
%       R-by-T matrix of graph spectral coefficients.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   The transform is computed as U'*X.
%
%   See also
%   GSPTB_iGFT



    Xhat = U'*X;



end
