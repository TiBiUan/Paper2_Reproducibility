function [RecError] = GSPTB_ComputeLaplacianReconstructionError(U1,L2)

    % L2 expressed in the U1 basis
    L2_1 = U1'*L2*U1;

    % Reconstructed Laplacian in the U1 basis
    L2_rec = U1*diag(diag(L2_1))*U1';

    % Error between reconstructed and actual Laplacians
    RecError = norm(L2-L2_rec,'fro')/norm(L2,'fro');
end