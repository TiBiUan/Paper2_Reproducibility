function [MeanIntensity] = GSPTB_ComputeMeanIntensity(ActiveTags,Weights)
%GSPTB_COMPUTEMEANINTENSITY Compute mean absolute expression intensity over active time points.
%
%   MeanIntensity = GSPTB_ComputeMeanIntensity(ActiveTags,Weights)
%
%   Description
%   -----------
%   Compute mean absolute expression intensity over active time points.
%
%   Inputs
%   ------
%   ActiveTags
%       R-by-T logical matrix identifying active coefficients.
%   Weights
%       R-by-T matrix of spectral coefficients.
%
%   Outputs
%   -------
%   MeanIntensity
%       R-element vector of mean absolute coefficient magnitudes; NaN for
%       harmonics that are never active.
%
%   API status
%   ----------
%   Public.
%
%   See also
%   GSPTB_ExtractDynamicMetrics, GSPTB_ComputePower

    
    % AbsWeights contains the mean absolute intensity value, kept as NaN if
    % a harmonic is never expressed
    MeanIntensity = NaN(size(ActiveTags,1),1);

    for m = 1:size(ActiveTags,1)
        is_active = ActiveTags(m,:);
    
        if any(is_active)
            MeanIntensity(m) = mean(abs(Weights(m,is_active)));
        end
    end
end
