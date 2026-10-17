function [Stability_CH,MU_CH,STD_CH] = GSPTB_ComputeHarmonicStability(U_R)
%GSPTB_COMPUTEHARMONICSTABILITY Compute population-wise regional stability of aligned harmonics.
%
%   [Stability_CH,MU_CH,STD_CH] = GSPTB_ComputeHarmonicStability(U_R)
%
%   Description
%   -----------
%   Compute population-wise regional stability of aligned harmonics.
%
%   Inputs
%   ------
%   U_R
%       1-by-S cell array of aligned R-by-R harmonic bases.
%
%   Outputs
%   -------
%   Stability_CH
%       R-by-R matrix of regional mean-to-standard-deviation ratios.
%   MU_CH
%       R-by-R matrix of regional harmonic means across subjects.
%   STD_CH
%       R-by-R matrix of regional harmonic standard deviations across subjects.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Rows index regions and columns index harmonics.
%
%   See also
%   GSPTB_ComputeHarmonicStabilityNull, GSPTB_ComputeCohortHarmonicDescriptors


    % Number of brain regions
    R = size(U_R{1},1);

    % Number of subjects
    S = length(U_R);

    % We build the matrix of a given harmonic c across all S subjects (size
    % R x S), inside CH_OI
    for c = 1:R

        % Initializes
        CH_OI = NaN(R,S);

        % Fills
        for s = 1:S
            CH_OI(:,s) = U_R{s}(:,c);
        end
            
        % Gets mean and STD of signal across subjects
        MU_CH(:,c) = mean(CH_OI,2);
        STD_CH(:,c) = std(CH_OI,[],2);
    end
    
    % Stability is the ratio
    Stability_CH = MU_CH./STD_CH;
end
