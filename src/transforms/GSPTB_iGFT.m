function [X] = GSPTB_iGFT(U,Xhat)
%GSPTB_IGFT Reconstruct graph signals from graph Fourier coefficients.
%
%   X = GSPTB_iGFT(U,Xhat)
%
%   Description
%   -----------
%   Reconstruct graph signals from graph Fourier coefficients.
%
%   Inputs
%   ------
%   U
%       R-by-R graph spectral basis.
%   Xhat
%       R-by-T matrix of graph spectral coefficients.
%
%   Outputs
%   -------
%   X
%       R-by-T reconstructed graph-signal matrix.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   The inverse transform is computed as U*Xhat.
%
%   See also
%   GSPTB_GFT



    X = U*Xhat;



end
