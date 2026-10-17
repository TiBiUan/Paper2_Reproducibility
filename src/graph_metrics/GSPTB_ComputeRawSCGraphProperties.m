function [Strength_RawSC,CC_RawSC,Eff_RawSC,BTW_RawSC,EIC_RawSC,PC_RawSC] = GSPTB_ComputeRawSCGraphProperties(A)
%GSPTB_COMPUTERAWSCGRAPHPROPERTIES Compute nodal graph-theoretical properties from structural connectivity.
%
%   [Strength_RawSC,CC_RawSC,Eff_RawSC,BTW_RawSC,EIC_RawSC,PC_RawSC] = GSPTB_ComputeRawSCGraphProperties(A)
%
%   Description
%   -----------
%   Compute nodal graph-theoretical properties from structural connectivity.
%
%   Inputs
%   ------
%   A
%       R-by-R weighted undirected structural connectivity matrix.
%
%   Outputs
%   -------
%   Strength_RawSC
%       R-element nodal strength vector.
%   CC_RawSC
%       R-element weighted clustering-coefficient vector.
%   Eff_RawSC
%       R-element nodal efficiency vector.
%   BTW_RawSC
%       R-element weighted betweenness-centrality vector.
%   EIC_RawSC
%       R-element eigenvector-centrality vector.
%   PC_RawSC
%       R-element participation-coefficient vector.
%
%   API status
%   ----------
%   Public.
%
%   Dependencies
%   ------------
%   strengths_und
%   clustering_coef_wu
%   efficiency_wei
%   betweenness_wei
%   eigenvector_centrality_und
%   modularity_und
%   participation_coef


    Strength_RawSC = strengths_und(A);
    CC_RawSC = clustering_coef_wu(A);
    Eff_RawSC = efficiency_wei(A,2);
    BTW_RawSC = betweenness_wei(1./A);
    EIC_RawSC = eigenvector_centrality_und(A);
    
    [Ci,~] = modularity_und(A);
    PC_RawSC = participation_coef(A,Ci);
end
