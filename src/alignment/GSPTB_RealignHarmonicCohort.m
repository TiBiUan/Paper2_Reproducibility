function [U_R,Lambda_R,idx_hung] = GSPTB_RealignHarmonicCohort(U,Lambda,s_ref)
%GSPTB_REALIGNHARMONICCOHORT Realign a cohort of subject-specific harmonic bases to a common reference.
%
%   [U_R,Lambda_R,idx_hung] = GSPTB_RealignHarmonicCohort(U,Lambda,s_ref)
%
%   Description
%   -----------
%   Realign a cohort of subject-specific harmonic bases to a common reference.
%
%   Inputs
%   ------
%   U
%       1-by-S cell array; U{s} is an R-by-R harmonic basis.
%   Lambda
%       1-by-S cell array; Lambda{s} contains R eigenvalues.
%   s_ref
%       Scalar index of the reference subject.
%
%   Outputs
%   -------
%   U_R
%       1-by-S cell array of reordered and sign-aligned harmonic bases.
%   Lambda_R
%       1-by-S cell array of eigenvalues reordered consistently with U_R.
%   idx_hung
%       1-by-S cell array of Hungarian assignment indices.
%
%   API status
%   ----------
%   Public.
%
%   Dependencies
%   ------------
%   GSPTB_RealignHarmonics_Hungarian
%
%   See also
%   GSPTB_RealignHarmonics_Hungarian, GSPTB_ComputeSortedHarmonics

    
    S = length(U);

    for s = 1:S
        [U_R{s},Lambda_R{s},idx_hung{s}] = GSPTB_RealignHarmonics_Hungarian(U{s},U{s_ref},Lambda{s});
    end
end
