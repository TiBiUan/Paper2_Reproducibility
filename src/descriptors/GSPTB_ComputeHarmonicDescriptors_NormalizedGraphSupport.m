function [NormalizedGraphSupport] = GSPTB_ComputeHarmonicDescriptors_NormalizedGraphSupport(u,A)
%GSPTB_COMPUTEHARMONICDESCRIPTORS_NORMALIZEDGRAPHSUPPORT Compute normalized graph support for one harmonic.
%
%   NormalizedGraphSupport = GSPTB_ComputeHarmonicDescriptors_NormalizedGraphSupport(u,A)
%
%   Description
%   -----------
%   Compute normalized graph support for one harmonic.
%
%   Inputs
%   ------
%   u
%       R-element harmonic vector.
%   A
%       R-by-R weighted structural connectivity matrix.
%
%   Outputs
%   -------
%   NormalizedGraphSupport
%       Scalar area under the normalized magnitude-versus-graph-distance profile.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Graph distances are computed from inverse edge weights.
%   Distances and absolute harmonic magnitudes are normalized by their maxima
%   before numerical integration.
%   If multiple regions share the maximum absolute harmonic value, the first
%   is used as the reference peak.
%
%   Dependencies
%   ------------
%   distance_wei
%
%   See also
%   GSPTB_ComputeHarmonicDescriptors


    % Distance on the graph from the structural connectivity matrix
    graph_dist = distance_wei(1./A);

    % Finds the maximum of the harmonic
    max_idx = find(abs(u) == max(abs(u)));

    % If there are two equal maxima, takes the first one
    max_idx = max_idx(1);

    % Initializes
    wavedist = NaN(length(u),1);

    % Looking at all regions
    for r = 1:length(u)

        % If we consider a region that is not the maximum signal one, then
        % we get the associated distance value
        if r ~= max_idx
            wavedist(r) = graph_dist(max_idx,r);
            wavedist(isinf(wavedist)) = NaN;

        % The distance is trivially 0 for our maximum
        else
            wavedist(r) = 0;
        end
    end

    % Normalizes the distances to min = 0, max = 1
    WaveletDistances = wavedist/max(wavedist);

    % Normalizes the magnitudes to min = 0, max = 1
    WaveletMagnitude = abs(u)/max(abs(u));

    % Sorts the distances in ascending order
    [~,id] = sort(WaveletDistances,'ascend');

    % Computes the area under the curve on the reordered distance and
    % magnitude data for trapz to function properly
    NormalizedGraphSupport = trapz(WaveletDistances(id),WaveletMagnitude(id));
end
