function [Sparsity,sSC,Localization,selvec] = GSPTB_ComputeHarmonicDescriptors_Sparsity_sSC_Localization(u,A,Dis,T)
%GSPTB_COMPUTEHARMONICDESCRIPTORS_SPARSITY_SSC_LOCALIZATION Compute sparsity, spanned connectivity, and anatomical localization for one harmonic.
%
%   [Sparsity,sSC,Localization,selvec] = GSPTB_ComputeHarmonicDescriptors_Sparsity_sSC_Localization(u,A,Dis,T)
%
%   Description
%   -----------
%   Compute sparsity, spanned connectivity, and anatomical localization for
%   one harmonic.
%
%   Inputs
%   ------
%   u
%       R-element harmonic vector.
%   A
%       R-by-R structural connectivity matrix.
%   Dis
%       R-by-R inter-regional Euclidean distance matrix.
%   T
%       Non-negative absolute-amplitude threshold.
%
%   Outputs
%   -------
%   Sparsity
%       Fraction of entries set to zero after thresholding.
%   sSC
%       Mean structural connectivity among pairs of retained regions.
%   Localization
%       Mean inverse Euclidean distance among pairs of retained regions.
%   selvec
%       R-element logical vector identifying retained regions.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Only upper-triangular region pairs are included in sSC and Localization.
%
%   Dependencies
%   ------------
%   jUpperTriMatToVec
%
%   See also
%   GSPTB_ComputeHarmonicDescriptors_Sparsity,
%   GSPTB_ComputeHarmonicDescriptors


    % We make the harmonic sparse according to threshold T
    u_sparse = u;
    u_sparse(abs(u)<T) = 0;

    SC_SelMat = logical(u_sparse*u_sparse');

    % Keeps track of which regions survived thresholding for the considered
    % harmonic and subject
    selvec = logical(u_sparse);

    % Sparsity
    Sparsity = sum(u_sparse==0)/length(u);

    % Spanned structural connectivity
    A_vec = jUpperTriMatToVec(A);
    sSC = mean(A_vec(jUpperTriMatToVec(SC_SelMat)));

    % Localization
    Dis_vec = jUpperTriMatToVec(Dis);
    Localization = mean(1./Dis_vec(jUpperTriMatToVec(SC_SelMat)));
end
