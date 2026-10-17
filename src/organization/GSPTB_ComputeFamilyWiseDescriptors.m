function [FamilyDescriptors] = GSPTB_ComputeFamilyWiseDescriptors(idx_fam,...
    Lambda,HarmDescriptors)
%GSPTB_COMPUTEFAMILYWISEDESCRIPTORS Aggregate harmonic descriptors within spectral families.
%
%   FamilyDescriptors = GSPTB_ComputeFamilyWiseDescriptors(idx_fam,Lambda,HarmDescriptors)
%
%   Description
%   -----------
%   Aggregate harmonic descriptors within spectral families.
%
%   Inputs
%   ------
%   idx_fam
%       R-element family-assignment vector.
%   Lambda
%       1-by-S cell array of R-element eigenvalue vectors.
%   HarmDescriptors
%       Structure of harmonic-wise descriptors.
%
%   Outputs
%   -------
%   FamilyDescriptors
%       Structure containing included harmonics, spectral span, size, mean index,
%       mean frequency, and family-averaged descriptors.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Families are assumed to be indexed consecutively from 1 to max(idx_fam).
%
%   See also
%   GSPTB_ExtractFamilies, GSPTB_ComputeCohortHarmonicDescriptors


    % Number of subjects
    S = length(Lambda);

    % Number of families
    F = max(idx_fam);

    % Number of harmonics
    R = length(idx_fam);

    % Creates spectral index
    spec_ind = 1:R;

    % We compute all the family-wise properties
    for f = 1:F
    
        % Included harmonic indices in each family
        FamilyDescriptors.IncludedHarmonics{f} = spec_ind(idx_fam == f);
    
        % Spectral span
        FamilyDescriptors.SpectralSpan(f) = max(spec_ind(idx_fam == f))-min(spec_ind(idx_fam == f));
    
        % Family size
        FamilyDescriptors.Size(f) = length(spec_ind(idx_fam == f));
    
        % Mean spectral index
        FamilyDescriptors.MeanSpectralIndex(f) = mean(spec_ind(idx_fam == f));

        % Mean spectral frequency
        for s = 1:S
            FamilyDescriptors.MeanSpectralFrequency(f,s) = mean(Lambda{s}(idx_fam == f));
        end

        % Absolute stability (F x n_regions)
        FamilyDescriptors.AbsoluteStability(f,:) = mean(abs(HarmDescriptors.Stability(:,idx_fam == f)),2);
    
        % Normalized graph support (F x S)
        FamilyDescriptors.GraphSupport(f,:) = mean(HarmDescriptors.GraphSupport(:,idx_fam == f),2);
        
        % Sparsity (F x S)
        FamilyDescriptors.Sparsity(f,:) = mean(HarmDescriptors.Sparsity(:,idx_fam == f),2);
    
        % Localization (F x S)
        FamilyDescriptors.Localization(f,:) = mean(HarmDescriptors.Localization(:,idx_fam == f),2);
    
        % Spanned SC (F x S)
        FamilyDescriptors.SpannedConnectivity(f,:) = mean(HarmDescriptors.SpannedConnectivity(:,idx_fam == f),2);
    end
end
