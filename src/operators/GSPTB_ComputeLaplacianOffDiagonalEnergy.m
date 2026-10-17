function [OffDiagEnergy] = GSPTB_ComputeLaplacianOffDiagonalEnergy(U1,L2)

    % L2 expressed in the U1 basis
    L2_1 = U1'*L2*U1;

    % We take the off-diagonal elements as a measure of "mismatch"
    OffDiagEnergy = norm((L2_1-diag(diag(L2_1))),'fro')^2/norm(L2_1,'fro')^2;
end