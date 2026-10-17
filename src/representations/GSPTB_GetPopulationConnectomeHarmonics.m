function [H] = GSPTB_GetPopulationConnectomeHarmonics(A)

    % Computations are carried for the input adjacency 
    L = GSPTB_ComputeNormalizedLaplacian(A);
    [U,Lambda] = GSPTB_ComputeSortedHarmonics(L);

    H.L = L;
    H.U = U;
    H.Lambda = Lambda;
end