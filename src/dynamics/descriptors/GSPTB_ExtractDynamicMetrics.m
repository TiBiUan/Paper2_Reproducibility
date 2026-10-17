function [EntryRate,Duration] = GSPTB_ExtractDynamicMetrics(BinWeights,TR)
%GSPTB_EXTRACTDYNAMICMETRICS Compute harmonic entry rate and mean active-period duration.
%
%   [EntryRate,Duration] = GSPTB_ExtractDynamicMetrics(BinWeights,TR)
%
%   Description
%   -----------
%   Compute harmonic entry rate and mean active-period duration.
%
%   Inputs
%   ------
%   BinWeights
%       R-by-T binary matrix of harmonic activity states.
%   TR
%       Sampling interval in seconds.
%
%   Outputs
%   -------
%   EntryRate
%       R-element vector of inactive-to-active transitions per second.
%   Duration
%       R-element vector of mean active-period durations in seconds.
%
%   API status
%   ----------
%   Public.
%
%   Dependencies
%   ------------
%   CAP_ComputeMetrics
%
%   See also
%   GSPTB_ComputeMeanIntensity, GSPTB_ComputePower


    % Number of frames to consider
    n_frames = size(BinWeights,2);

    [~,~,Number,Avg_Duration] = CAP_ComputeMetrics(single(BinWeights),1,TR,n_frames);

    % Entry rate in [1/s]
    EntryRate = Number(:,3)/n_frames/TR;

    % Duration in [s]
    Duration = Avg_Duration(:,3);
end
