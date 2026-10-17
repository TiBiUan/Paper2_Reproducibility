function [R2] = GSPTB_ComputeR2_Covariate(y_fit,y,cova)
%GSPTB_COMPUTER2_COVARIATE Compute covariate-adjusted coefficient of determination.
%
%   R2 = GSPTB_ComputeR2_Covariate(y_fit,y,cova)
%
%   Description
%   -----------
%   Compute covariate-adjusted coefficient of determination.
%
%   Inputs
%   ------
%   y_fit
%       Vector of fitted values.
%   y
%       Vector of observed values.
%   cova
%       N-by-P matrix of covariates of no interest.
%
%   Outputs
%   -------
%   R2
%       Coefficient of determination computed between covariate-residualized
%       fitted and observed values.
%
%   API status
%   ----------
%   Public.
%
%   Dependencies
%   ------------
%   y_regress_ss
%
%   See also
%   GSPTB_ComputeR2


    [~,res_fit] = y_regress_ss(y_fit,cova);
    [~,res] = y_regress_ss(y,cova);

    R2 = 1 - sum((res_fit-res).^2)/sum((mean(res)-res).^2);
end
