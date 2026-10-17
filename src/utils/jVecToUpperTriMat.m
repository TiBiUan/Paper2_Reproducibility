function m=jVecToUpperTriMat(v,size_m)
%JVECTOUPPERTRIMAT Insert a vector into the strict upper triangle of a square matrix.
%
%   m = jVecToUpperTriMat(v,size_m)
%
%   Description
%   -----------
%   Insert a vector into the strict upper triangle of a square matrix.
%
%   Inputs
%   ------
%   v
%       Vector of upper-triangular values.
%   size_m
%       Scalar side length of the output square matrix.
%
%   Outputs
%   -------
%   m
%       size_m-by-size_m matrix with v in the strict upper triangle and zeros
%       elsewhere.
%
%   API status
%   ----------
%   External.
%
%   Notes
%   -----
%   External laboratory utility authored by Dimitri Van De Ville (2015).
%   The lower triangle is not filled automatically.
%
%   See also
%   jUpperTriMatToVec

% get indices of upper triangular part (Peter Acklam's trick)
idx = find(triu(ones(size_m), 1));

m=zeros(size_m,size_m);
m(idx)=v(:);
%m=m+m.';
