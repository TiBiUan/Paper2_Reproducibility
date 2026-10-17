function [d_P] = GSPTB_ComputeSubspaceProjectionMatrixDistance(U1,U2)

    if size(U1,1) ~= size(U2,1) || size(U1,2) ~= size(U2,2)
        error('U1 and U2 must define subspaces of equal dimensions!');
    end

    % Dimension of subspace
    k = size(U1,2);

    % Number of harmonics in each subspace
    P1 = U1*U1';
    P2 = U2*U2';
    
    % Frobenius norm of the difference between both projection matrices
    d_P = norm(P1-P2,'fro')/sqrt(2*k);
end