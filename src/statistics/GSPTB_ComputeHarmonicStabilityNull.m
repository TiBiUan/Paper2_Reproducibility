function [Stability_SignElements] = GSPTB_ComputeHarmonicStabilityNull(Stability_CH,U_R,n_null_stability,alpha_stab,n_tests_stab)
%GSPTB_COMPUTEHARMONICSTABILITYNULL Identify region-harmonic stability values exceeding a permutation-derived null threshold.
%
%   Stability_SignElements = GSPTB_ComputeHarmonicStabilityNull(Stability_CH,U_R,n_null_stability,alpha_stab,n_tests_stab)
%
%   Description
%   -----------
%   Identify region-harmonic stability values exceeding a permutation-derived
%   null threshold.
%
%   Inputs
%   ------
%   Stability_CH
%       R-by-R observed regional stability matrix.
%   U_R
%       1-by-S cell array of aligned R-by-R harmonic bases.
%   n_null_stability
%       Number of null realizations.
%   alpha_stab
%       Family-wise significance level expressed in percent units for prctile.
%   n_tests_stab
%       Number of parallel tests used for Bonferroni correction.
%
%   Outputs
%   -------
%   Stability_SignElements
%       R-by-R logical matrix marking observed stability values above the
%       corrected null threshold.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Within each null realization, regional entries are independently permuted
%   for every subject and harmonic.
%
%   Dependencies
%   ------------
%   GSPTB_ComputeHarmonicStability
%
%   See also
%   GSPTB_ComputeHarmonicStability


    % Number of brain regions
    R = size(Stability_CH,1);

    % Number of subjects
    S = length(U_R);
    
    % For null stability, we will concatenate the outputs across regions for each
    % harmonic
    Stability_Null = NaN(R,R,n_null_stability);
    
    % We randomize the regions
    for n = 1:n_null_stability
    
        n
    
        % We will randomize U_R_NullStab differently every time
        U_R_NullStab = U_R;
    
        % I randomize each harmonic from each subject differently
        for s = 1:S
            for c = 1:R
                
                % Random ordering
                idx_rand = randperm(R);
                U_R_NullStab{s}(:,c) = U_R_NullStab{s}(idx_rand,c);
            end
        end
    
        % Computes null stability
        [Stability_Null(:,:,n)] = abs(GSPTB_ComputeHarmonicStability(U_R_NullStab)');
    end
    
    % Gets the result when concatenating across regions
    % Resulting size: n_harmonics x pooled_null_observations
    Stability_Null = reshape(Stability_Null,R,R*n_null_stability);
    
    % We can compute a significance threshold and extract the set of
    % significant elements (Supplementary Figure 2D)
    T_Stability = prctile(Stability_Null,100-alpha_stab/n_tests_stab,2);
    Stability_SignElements = (abs(Stability_CH) > repmat(T_Stability,1,R)');
end
