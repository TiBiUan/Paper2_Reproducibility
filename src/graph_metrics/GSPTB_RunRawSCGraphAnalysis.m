function [RawSCAnalysis] = GSPTB_RunRawSCGraphAnalysis(Adj,SpecPart)
    
    % Number of subjects
    S = length(Adj);

    % Number of regions
    R = size(Adj{1},1);

    % Computes graph theoretical properties on raw adjacency matrices
    Strength_RawSC = NaN(R,S);
    CC_RawSC = NaN(R,S);
    Eff_RawSC = NaN(R,S);
    BTW_RawSC = NaN(R,S);
    EIC_RawSC = NaN(R,S);
    PC_RawSC = NaN(R,S);

    for s = 1:S
    
        s
    
        [Strength_RawSC(:,s),CC_RawSC(:,s),Eff_RawSC(:,s),BTW_RawSC(:,s),...
            EIC_RawSC(:,s),PC_RawSC(:,s)] = GSPTB_ComputeRawSCGraphProperties(Adj{s});
    end

    RawSCAnalysis.Strength = Strength_RawSC;
    RawSCAnalysis.ClusteringCoefficient = CC_RawSC;
    RawSCAnalysis.Efficiency = Eff_RawSC;
    RawSCAnalysis.BetweennessCentrality = BTW_RawSC;
    RawSCAnalysis.EigenvectorCentrality = EIC_RawSC;
    RawSCAnalysis.ParticipationCoefficient = PC_RawSC;

    % Correlations with spectral participation
    % Size S x 6, with the graph-theoretical properties ordered as
    % Strength | CC | Efficiency | Betweenness | Eigen centrality | Part coef
    RawSCAnalysis.CorrelationWithSpectralParticipation = ...
        [diag(corr(SpecPart,Strength_RawSC)),...
        diag(corr(SpecPart,CC_RawSC)),...
        diag(corr(SpecPart,Eff_RawSC)),...
        diag(corr(SpecPart,EIC_RawSC)),...
        diag(corr(SpecPart,BTW_RawSC)),...
        diag(corr(SpecPart,PC_RawSC))];
end