function [HarmonicDescriptors] = ...
    GSPTB_ComputeCohortHarmonicDescriptors(U,Lambda,A,Deg,T_Sparsity,Dis,is_normalize,gamma)
%GSPTB_COMPUTECOHORTHARMONICDESCRIPTORS Compute harmonic descriptors for every subject and harmonic in a cohort.
%
%   HarmonicDescriptors = GSPTB_ComputeCohortHarmonicDescriptors(U,Lambda,A,Deg,T_Sparsity,Dis,is_normalize,gamma)
%
%   Description
%   -----------
%   Compute harmonic descriptors for every subject and harmonic in a cohort.
%
%   Inputs
%   ------
%   U
%       1-by-S cell array of R-by-R harmonic bases.
%   Lambda
%       1-by-S cell array of R-element eigenvalue vectors.
%   A
%       1-by-S cell array of R-by-R structural connectivity matrices.
%   Deg
%       1-by-S cell array of R-by-R degree matrices.
%   T_Sparsity
%       Scalar threshold used to define supra-threshold harmonic signal.
%   Dis
%       R-by-R matrix of inter-regional Euclidean distances.
%   is_normalize
%       Logical flag controlling Frobenius normalization of reconstructed
%       connectivity.
%   gamma
%       Scalar resolution parameter for modularity estimation.
%
%   Outputs
%   -------
%   HarmonicDescriptors
%       Structure containing stability, graph support, sparsity, spanned
%       connectivity, localization, zero crossings, reconstructed connectivity,
%       efficiency, modularity, selection vectors, and spectral participation.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Population stability is computed once across the aligned cohort;
%   subject-wise descriptors are then evaluated for every harmonic.
%
%   Dependencies
%   ------------
%   GSPTB_ComputeHarmonicStability
%   GSPTB_ComputeHarmonicDescriptors
%
%   See also
%   GSPTB_ComputeHarmonicDescriptors, GSPTB_ComputeFamilyWiseDescriptors


    % Number of subjects
    S = length(U);

    % Number of regions
    R = size(U{1},1);

    % Computation of harmonic stability
    [HarmonicDescriptors.Stability,HarmonicDescriptors.StabNum,...
        HarmonicDescriptors.StabDenom] = GSPTB_ComputeHarmonicStability(U);
    
    % Same process for each subject and harmonic...
    for s = 1:S

        s
    
        for c = 1:R
    
            c
        
            % Computation of all descriptors
            [HarmonicDescriptors.GraphSupport(s,c),HarmonicDescriptors.Sparsity(s,c),...
                HarmonicDescriptors.SpannedConnectivity(s,c),HarmonicDescriptors.Localization(s,c),...
                HarmonicDescriptors.ZeroCrossings(s,c),HarmonicDescriptors.ApproximatedSC{s}{c},...
                HarmonicDescriptors.Efficiency(s,c),HarmonicDescriptors.Modularity(s,c),...
                HarmonicDescriptors.SelectionVectors(:,s,c)] = ...
                GSPTB_ComputeHarmonicDescriptors(U{s}(:,c),Lambda{s}(c),...
                A{s},Deg{s},T_Sparsity,Dis,is_normalize,gamma);
        end
    end

    HarmonicDescriptors.SpectralParticipation = nansum(HarmonicDescriptors.SelectionVectors,3);
end
