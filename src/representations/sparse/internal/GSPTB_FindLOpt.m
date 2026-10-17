function [L_opt,RMSE_train] = GSPTB_FindLOpt(Dico,A_train)
%GSPTB_FINDLOPT Select the orthogonal-matching-pursuit sparsity level from an RMSE knee point.
%
%   [L_opt,RMSE_train] = GSPTB_FindLOpt(Dico,A_train)
%
%   Description
%   -----------
%   Select the orthogonal-matching-pursuit sparsity level from an RMSE knee
%   point.
%
%   Inputs
%   ------
%   Dico
%       N-by-K dictionary matrix.
%   A_train
%       N-by-T training signal matrix.
%
%   Outputs
%   -------
%   L_opt
%       Selected maximum number of active atoms.
%   RMSE_train
%       K-element reconstruction-RMSE curve for candidate sparsity levels 1:K.
%
%   API status
%   ----------
%   Internal.
%
%   Notes
%   -----
%   Internal helper for sparse representation estimation.
%
%   Dependencies
%   ------------
%   mexOMP
%   knee_pt
%
%   See also
%   GSPTB_RunSparseActivityPrediction


    % Number of atoms
    n_atoms = size(Dico,2);

    % Number of time points
    T = size(A_train,2);

    % Size of the vector to predict at each time point
    N = size(A_train,1);

    % Candidate values to examine
    L = 1:n_atoms;

    % Computes the weights for each candidate "sparsity" level
    for l = 1:length(L)

        l

        % The parameter L must be part of the "param" struture for use with
        % the dictionary learning (SPAMS) toolbox
        param.L = L(l);

        % Gets the weights for the given L
        Weights{l} = full(mexOMP(A_train,single(Dico),param));

        % Computes the error
        RMSE_train(l) = sqrt(sum(sum(((A_train-Dico*Weights{l}).^2)))/(N*T));
    end

    % We wish to determine where, in the error curve, we have no more
    % "improvement"
    knee = knee_pt(RMSE_train);

    L_opt = L(knee);
end
