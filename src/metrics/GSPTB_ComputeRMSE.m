function [RMSE] = GSPTB_ComputeRMSE(y_fit,y)
%GSPTB_COMPUTERMSE Compute root mean square error.
%
%   RMSE = GSPTB_ComputeRMSE(y_fit,y)
%
%   Description
%   -----------
%   Compute root mean square error.
%
%   Inputs
%   ------
%   y_fit
%       Vector of fitted values.
%   y
%       Observed vector with the same number of elements.
%
%   Outputs
%   -------
%   RMSE
%       Root mean square error.
%
%   API status
%   ----------
%   Public.
%
%   See also
%   GSPTB_ComputeR2


    N = length(y);

    RMSE = sqrt(sum((y_fit-y).^2)/N);
end
