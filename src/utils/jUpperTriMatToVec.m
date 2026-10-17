function v=jUpperTriMatToVec(m,varargin)
%JUPPERTRIMATTOVEC Convert the upper-triangular part of a matrix to a vector.
%
%   v = jUpperTriMatToVec(m,offset)
%
%   Description
%   -----------
%   Convert the upper-triangular part of a matrix to a vector.
%
%   Inputs
%   ------
%   m
%       Input matrix.
%   offset
%       Optional diagonal offset passed to triu; defaults to 1.
%
%   Outputs
%   -------
%   v
%       Vector containing selected upper-triangular matrix entries.
%
%   API status
%   ----------
%   External.
%
%   Notes
%   -----
%   External laboratory utility originally authored by Jonas Richiardi (2009)
%   and revised by Dimitri Van De Ville (2015).
%
%   See also
%   jVecToUpperTriMat

switch nargin
    case 1
        offset=1;
    case 2
        offset=varargin{1};
end

% get indices of upper triangular part (Peter Acklam's trick)
%[m_i m_j] = find(triu(ones(size(m)), offset));
idx = find(triu(ones(size(m)), offset));

v=m(idx);

% copy to vector
%v=zeros(numel(m_i),1);
%for v_idx=1:numel(m_i)
%    v(v_idx)=m(m_i(v_idx),m_j(v_idx));
%end
