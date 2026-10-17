function [U_new,Lambda_new,idx] = GSPTB_RealignHarmonics_Hungarian(U,U_ref,Lambda)
%GSPTB_REALIGNHARMONICS_HUNGARIAN Align one harmonic basis to a reference basis using the Hungarian algorithm.
%
%   [U_new,Lambda_new,idx] = GSPTB_RealignHarmonics_Hungarian(U,U_ref,Lambda)
%
%   Description
%   -----------
%   Align one harmonic basis to a reference basis using the Hungarian
%   algorithm.
%
%   Inputs
%   ------
%   U
%       R-by-R harmonic basis to align.
%   U_ref
%       R-by-R reference harmonic basis.
%   Lambda
%       R-element eigenvalue vector associated with U.
%
%   Outputs
%   -------
%   U_new
%       R-by-R reordered and sign-aligned harmonic basis.
%   Lambda_new
%       R-element eigenvalue vector reordered according to idx.
%   idx
%       R-element assignment vector returned by the Hungarian algorithm.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Assignment costs are 1 minus the absolute Pearson correlation between
%   harmonic spatial profiles.
%   After reordering, each harmonic sign is chosen to yield non-negative
%   correlation with its reference counterpart.
%
%   Dependencies
%   ------------
%   munkres
%
%   See also
%   GSPTB_RealignHarmonicCohort


    % Atlas size
    R = size(U,1);
    
    % Extracts the best ordering
    idx = munkres((1-abs(corr(U_ref,U))));

    % Reorders the harmonics
    U_new = U(:,idx);

    % Adjusts their signs
    sign_corr = sign(diag(corr(U_ref,U_new)));
    sign_corr(sign_corr == 0) = 1;
    U_new = U_new.*repmat(sign_corr',R,1);

    % Reorders the eigenvalues
    Lambda_new = Lambda(idx);
end
