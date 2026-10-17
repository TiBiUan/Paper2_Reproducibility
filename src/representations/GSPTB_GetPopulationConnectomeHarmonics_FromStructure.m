function [U,Lambda] = GSPTB_GetPopulationConnectomeHarmonics_FromStructure(A_struct)

    

    % Computations are carried for the input adjacency 
    L = GSPTB_ComputeNormalizedLaplacian(A);
    [U,Lambda] = GSPTB_ComputeSortedHarmonics(L);
end