function Power = GSPTB_ComputePower(Weights,FD,T_mot)
%GSPTB_COMPUTEPOWER Compute harmonic power over retained low-motion time points.
%
%   Power = GSPTB_ComputePower(Weights,FD,T_mot)
%
%   Description
%   -----------
%   Compute harmonic power over retained low-motion time points.
%
%   Inputs
%   ------
%   Weights
%       R-by-T matrix of sparse spectral coefficients.
%   FD
%       T-element framewise-displacement vector.
%   T_mot
%       Scalar framewise-displacement threshold.
%
%   Outputs
%   -------
%   Power
%       R-element vector of mean squared coefficients over retained time points.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Zero coefficients are included, so power jointly reflects expression
%   magnitude, duration, and recruitment frequency.
%   Time points are retained when FD is less than or equal to T_mot.
%
%   See also
%   GSPTB_ComputeMeanIntensity


    % Number of harmonics
    % R = size(Weights,1);

    FD = FD(:);

    if size(Weights,2) ~= numel(FD)
        error(['The number of columns in Weights must equal the ' ...
            'number of elements in FD!']);
    end

    % Is a given harmonic active?
    % is_harm_active = (Weights ~= 0);

    % Is movement sufficiently low?
    is_retained_TP = FD <= T_mot;

    % Are both conditions satisfied at once?
    % is_retained_dp = is_motion_ok;

    % Power = NaN(R,1);
    % 
    % for c = 1:R
    %     CurrentWeights = Weights(c,is_retained_dp(c,:));
    %     Power(c) = mean(abs(CurrentWeights).^2);
    % end

    % Mean squared coefficient across retained time points
    Power = mean(Weights(:,is_retained_TP).^2,2);
end
