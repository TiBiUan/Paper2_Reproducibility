function [Net_ZC] = GSPTB_ComputeHarmonicDescriptors_NetZC(u,A,T)
%GSPTB_COMPUTEHARMONICDESCRIPTORS_NETZC Compute the weighted net zero-crossing of one harmonic.
%
%   Net_ZC = GSPTB_ComputeHarmonicDescriptors_NetZC(u,A,T)
%
%   Description
%   -----------
%   Compute the weighted net zero-crossing of one harmonic.
%
%   Inputs
%   ------
%   u
%       R-element harmonic vector.
%   A
%       R-by-R structural connectivity matrix.
%   T
%       Scalar threshold applied to squared harmonic amplitudes.
%
%   Outputs
%   -------
%   Net_ZC
%       Sum of structural edge weights linking retained regions with opposite
%       harmonic signs.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Thresholding is intentionally applied to u.^2 to reproduce the formulation
%   of Sipes et al. (2024).
%   Only upper-triangular edge entries are counted.
%
%   Dependencies
%   ------------
%   jUpperTriMatToVec
%
%   See also
%   GSPTB_ComputeHarmonicDescriptors_Sparsity,
%   GSPTB_ComputeHarmonicDescriptors


    % Removes low values following the definition in Sipes et al. (2024)
    u_sparse = u;
    u_sparse(u.^2<=T)=0;

    SC_SelMat = (u_sparse*u_sparse')<0;

    % Computes net zero-crossing
    A_vec = jUpperTriMatToVec(A);
    Net_ZC = sum(A_vec(jUpperTriMatToVec(SC_SelMat)));

end
