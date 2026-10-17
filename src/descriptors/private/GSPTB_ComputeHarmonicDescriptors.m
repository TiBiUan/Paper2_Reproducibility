function [NormGraphSupport,Sparsity,sSC,Localization,Net_ZC,A_appr,Eff,Q,selvec] = ...
    GSPTB_ComputeHarmonicDescriptors(u,lam,A,Deg,T,Dis,is_normalize,gamma)
%GSPTB_COMPUTEHARMONICDESCRIPTORS Compute the complete descriptor set for one harmonic.
%
%   [NormGraphSupport,Sparsity,sSC,Localization,Net_ZC,A_appr,Eff,Q,selvec] = GSPTB_ComputeHarmonicDescriptors(u,lam,A,Deg,T,Dis,is_normalize,gamma)
%
%   Description
%   -----------
%   Compute the complete descriptor set for one harmonic.
%
%   Inputs
%   ------
%   u
%       R-element harmonic vector.
%   lam
%       Scalar eigenvalue associated with u.
%   A
%       R-by-R structural connectivity matrix.
%   Deg
%       R-by-R degree matrix.
%   T
%       Scalar harmonic threshold.
%   Dis
%       R-by-R inter-regional Euclidean distance matrix.
%   is_normalize
%       Logical flag controlling Frobenius normalization of A_appr.
%   gamma
%       Scalar modularity resolution parameter.
%
%   Outputs
%   -------
%   NormGraphSupport
%       Normalized graph-support measure.
%   Sparsity
%       Fraction of thresholded zero entries.
%   sSC
%       Mean structural connectivity among retained regions.
%   Localization
%       Mean inverse Euclidean distance among retained regions.
%   Net_ZC
%       Weighted net zero-crossing.
%   A_appr
%       Harmonic-wise reconstructed structural connectivity matrix.
%   Eff
%       Global efficiency of A_appr.
%   Q
%       Modularity of A_appr after adding unit self-connections.
%   selvec
%       Logical vector of retained regions.
%
%   API status
%   ----------
%   Internal.
%
%   Notes
%   -----
%   Internal orchestration helper used by the cohort-level descriptor routine.
%   Unit self-connections are added only for modularity to reproduce the
%   published implementation.
%
%   Dependencies
%   ------------
%   GSPTB_ComputeHarmonicDescriptors_NormalizedGraphSupport
%   GSPTB_ComputeHarmonicDescriptors_Sparsity_sSC_Localization
%   GSPTB_ComputeHarmonicDescriptors_NetZC
%   GSPTB_EstimateHarmonicSC
%   efficiency_wei
%   modularity_und
%
%   See also
%   GSPTB_ComputeCohortHarmonicDescriptors


    % Computes normalized graph support
    NormGraphSupport = GSPTB_ComputeHarmonicDescriptors_NormalizedGraphSupport(u,A);

    % Computes sparsity and spanned structural connectivity; also keeps
    % track of the regions surviving thresholding for spectral
    % participation assessments
    [Sparsity,sSC,Localization,selvec] = GSPTB_ComputeHarmonicDescriptors_Sparsity_sSC_Localization(u,A,Dis,T);

    % Computes net zero-crossing
    Net_ZC = GSPTB_ComputeHarmonicDescriptors_NetZC(u,A,T);
    
    % Generates the underlying equivalent structural connectivity matrix
    % and extracts graph properties
    A_appr = GSPTB_EstimateHarmonicSC(u,lam,Deg,is_normalize);

    % Following the original Paper implementation, unit self-connections are
    % added only for modularity computation. This convention is required to
    % reproduce the published modularity results
    A_appr_Q = A_appr;

    for r = 1:length(u)
        A_appr_Q(r,r) = 1;
    end

    % Efficiency and modularity
    Eff = efficiency_wei(A_appr);
    [~,Q] = modularity_und(A_appr_Q,gamma);
end
