function [Families] = GSPTB_ExtractFamilies(CrossSubjSimMat,FamOrg)
%GSPTB_EXTRACTFAMILIES Identify and order harmonic families by hierarchical clustering.
%
%   Families = GSPTB_ExtractFamilies(CrossSubjSimMat,FamOrg)
%
%   Description
%   -----------
%   Identify and order harmonic families by hierarchical clustering.
%
%   Inputs
%   ------
%   CrossSubjSimMat
%       R-by-R cross-subject harmonic similarity matrix.
%   FamOrg
%       Configuration structure defining ClusteringType, ClusteringDistance,
%       CuttingDistance, and CuttingCriterion.
%
%   Outputs
%   -------
%   Families
%       Structure containing family indices, number of families, reordered
%       harmonic order, and cumulative family boundaries.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Similarity profiles are z-scored before hierarchical clustering.
%   Families are relabeled by ascending mean harmonic index.
%
%   Dependencies
%   ------------
%   linkage
%   cluster
%   zscore
%
%   See also
%   GSPTB_ComputeFamilyDeltaSim, GSPTB_ComputeFamilyWiseDescriptors

    
    % Number of harmonics
    R = size(CrossSubjSimMat,1);

    % We run clustering on the z-scored version of the cross-subject similarity
    % Rows of the input matrix are the data points, hence the transpose
    Z = linkage(zscore(CrossSubjSimMat)',FamOrg.ClusteringType,FamOrg.ClusteringDistance);

    % Regional assignment indices for families, not yet reordered
    id_fams = cluster(Z,'cutoff',FamOrg.CuttingDistance,'criterion',FamOrg.CuttingCriterion);
    
    % Number of families at hand
    n_families = length(unique(id_fams));
    
    % Spectral indices in linear progression
    spec_ind = 1:R;
    
    % Computes mean spectral index across families
    mean_specid = NaN(n_families,1);

    for f = 1:n_families
        mean_specid(f) = mean(spec_ind(id_fams == f));
    end
    
    % Reorders the families as a function of mean spectral index
    [~,ind_fam_msi] = sort(mean_specid,'ascend');
    
    % We want to create a new index vector that tags harmonics to a family, so
    % that index 1 is the lowest spectral index
    id_fams_reordered = NaN(R,1);
    
    for f = 1:n_families
        id_fams_reordered(id_fams == ind_fam_msi(f)) = f;
    end

    perm = zeros(R,1);
    k = 1;
    
    for f = 1:n_families
        idx = sort(find(id_fams_reordered == f));
        perm(k:k+numel(idx)-1) = idx;
        k = k + numel(idx);
    end

    Families.Indices = id_fams_reordered;
    Families.NumFamilies = n_families;
    Families.Order = perm;
    Families.Boundaries = cumsum(histcounts(id_fams_reordered,1:n_families+1));
end
