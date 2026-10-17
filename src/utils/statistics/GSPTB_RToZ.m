function [z] = GSPTB_RToZ(R)
%GSPTB_RTOZ Apply the Fisher r-to-z transformation.
%
%   z = GSPTB_RToZ(R)
%
%   Description
%   -----------
%   Apply the Fisher r-to-z transformation.
%
%   Inputs
%   ------
%   R
%       Scalar, vector, or matrix of correlation coefficients.
%
%   Outputs
%   -------
%   z
%       Element-wise inverse hyperbolic tangent of R.
%
%   API status
%   ----------
%   Public.
%
%   Notes
%   -----
%   Values equal to plus or minus one map to infinite z scores.


    z = atanh(R);

end
