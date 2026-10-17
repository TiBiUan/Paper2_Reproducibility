function [SurrogateSparsityResults] = GSPTB_RunSurrogateSparsityAnalysis(Adj,T_Sparsity,s_ref)
%GSPTB_RUNSURROGATESPARSITYANALYSIS Assess harmonic sparsity under altered structural-connectivity constructions.
%
%   SurrogateSparsityResults = GSPTB_RunSurrogateSparsityAnalysis(Adj,T_Sparsity,s_ref)
%
%   Description
%   -----------
%   Assess harmonic sparsity under altered structural-connectivity
%   constructions.
%
%   Inputs
%   ------
%   Adj
%       1-by-S cell array of R-by-R structural connectivity matrices.
%   T_Sparsity
%       Scalar absolute-amplitude threshold.
%   s_ref
%       Reference-subject index for surrogate harmonic alignment.
%
%   Outputs
%   -------
%   SurrogateSparsityResults
%       Structure containing surrogate sparsity matrices, unaligned surrogate
%       bases, surrogate labels, and surrogate Laplacians.
%
%   API status
%   ----------
%   Public.
%
%   Dependencies
%   ------------
%   GSPTB_CreateSurrogateSCHarmonics
%   GSPTB_RealignHarmonics_Hungarian
%   GSPTB_ComputeHarmonicDescriptors_Sparsity
%
%   See also
%   GSPTB_CreateSurrogateSCHarmonics


    % Number of subjects
    S = length(Adj);

    % Number of regions
    R = size(Adj{1},1);

    % Surrogate harmonics for the reference subject
    [~,U_Surr_Ref,~] = GSPTB_CreateSurrogateSCHarmonics(Adj{s_ref});
    
    Ln_Surr = cell(1,S);
    U_Surr = cell(1,S);
    U_Surr_R = cell(1,S);

    Sparsity_Surr = cell(1,length(U_Surr_Ref));

    % Creates altered versions of the matrices and computes sparsity
    for s = 1:S
    
        s
    
        [Ln_Surr{s},U_Surr{s},Surr_Types] = GSPTB_CreateSurrogateSCHarmonics(Adj{s});

        U_Surr_R{s} = cell(1,length(U_Surr{s}));

        for i = 1:length(U_Surr{s})
            U_Surr_R{s}{i} = GSPTB_RealignHarmonics_Hungarian(U_Surr{s}{i},U_Surr_Ref{i},zeros(R,1));
            
            for c = 1:R
                Sparsity_Surr{i}(c,s) = GSPTB_ComputeHarmonicDescriptors_Sparsity(U_Surr_R{s}{i}(:,c),T_Sparsity);
            end
        end
    end

    SurrogateSparsityResults.Sparsity_Surr = Sparsity_Surr;
    SurrogateSparsityResults.U_Surr = U_Surr;
    SurrogateSparsityResults.Surr_Types = Surr_Types;
    SurrogateSparsityResults.Ln_Surr = Ln_Surr;
end
