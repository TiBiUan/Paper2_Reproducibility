function [R2] = GSPTB_ComputeR2(y_fit,y_actual)
%GSPTB_COMPUTER2 Compute the coefficient of determination between fitted and observed values.
%
%   R2 = GSPTB_ComputeR2(y_fit,y_actual)
%
%   Description
%   -----------
%   Compute the coefficient of determination between fitted and observed
%   values.
%
%   Inputs
%   ------
%   y_fit
%       Vector of fitted values.
%   y_actual
%       Observed vector with the same number of elements as y_fit.
%
%   Outputs
%   -------
%   R2
%       Coefficient of determination, 1 minus residual sum of squares divided by
%       total sum of squares.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   R2 may be negative for out-of-sample predictions.
%
%   See also
%   GSPTB_ComputeR2_Covariate, GSPTB_ComputeRMSE

    R2 = 1 - sum((y_fit-y_actual).^2)/sum((mean(y_actual)-y_actual).^2);
end
