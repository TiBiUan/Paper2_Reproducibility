function [A_Pop] = GSPTB_GetPopulationAdjacencyMatrices(Data,Cfg,id)

    % Computes the two population weight estimates
    A_All_All = GSPTB_ComputePopulationAverage(Data.Adj,id);
    A_All_NZ = GSPTB_ComputePopulationAverageNonZero(Data.Adj,id);

    % Computes the consensus mask and corresponding
    % consensus-edge population matrices
    G = GSPTB_ComputePopulationConsensusMask(Data.Adj,...
        id,Data.RegionDistances,Data.HemiID,Cfg.Population.NumBins);

    A_Cons_All = G.*A_All_All;
    A_Cons_NZ = G.*A_All_NZ;

    % Creates a structure to return that contains all the subcases
    A_Pop.AllEdges.AllValues = A_All_All;
    A_Pop.AllEdges.NZValues = A_All_NZ;

    A_Pop.ConsensusEdges.AllValues = A_Cons_All;
    A_Pop.ConsensusEdges.NZValues = A_Cons_NZ;
end