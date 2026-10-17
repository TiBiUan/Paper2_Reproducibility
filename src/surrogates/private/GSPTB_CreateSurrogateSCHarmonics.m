function [Ln,U,Surrogate_Types] = GSPTB_CreateSurrogateSCHarmonics(A)
%GSPTB_CREATESURROGATESCHARMONICS Generate structural-connectivity surrogates and their harmonic bases.
%
%   [Ln,U,Surrogate_Types] = GSPTB_CreateSurrogateSCHarmonics(A)
%
%   Description
%   -----------
%   Generate structural-connectivity surrogates and their harmonic bases.
%
%   Inputs
%   ------
%   A
%       R-by-R weighted undirected structural connectivity matrix.
%
%   Outputs
%   -------
%   Ln
%       1-by-4 cell array of normalized Laplacians.
%   U
%       1-by-4 cell array of sorted harmonic bases.
%   Surrogate_Types
%       1-by-4 cell array of descriptive surrogate labels.
%
%   API status
%   ----------
%   Internal.
%
%   Notes
%   -----
%   Internal helper generating original, binary, degree-preserving randomized,
%   and binary randomized conditions.
%
%   Dependencies
%   ------------
%   randmio_und_connected
%   GSPTB_ComputeNormalizedLaplacian
%   GSPTB_ComputeSortedHarmonics
%
%   See also
%   GSPTB_RunSurrogateSparsityAnalysis


    % Binary matrix
    A_binary = double(logical(A));
    
    % Randomized matrix
    A_rand = randmio_und_connected(A,20);
    
    % Binary and randomized matrix
    A_rand_binary = double(logical(A_rand));

    % Contains all sets of adjacency matrices
    A_set = {A,A_binary,A_rand,A_rand_binary};
    Surrogate_Types = {'Normal','Binary','Degree-preserving randomized','Binary & Degree-preserving randomized'};

    % Pre-allocates
    Ln = cell(1,length(A_set));
    U = cell(1,length(A_set));

    for i = 1:length(A_set)
        Ln{i} = GSPTB_ComputeNormalizedLaplacian(A_set{i});
        U{i} = GSPTB_ComputeSortedHarmonics(Ln{i});
    end
end
