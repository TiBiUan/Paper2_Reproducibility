function [SubSim,Angles,SV] = GSPTB_ComputeSubspacePrincipalAngles(U1,U2)

    if size(U1,1) ~= size(U2,1) || size(U1,2) ~= size(U2,2)
        error('U1 and U2 must define subspaces of equal dimensions!');
    end

    % Number of harmonics in each subspace
    n_sub = size(U1,2);

    % Singular values in a vector
    SV = svd(U1'*U2);

    % Angles matched to the singular values
    Angles = acos(min(max(SV,-1),1));
    
    % Subspace similarity: at most equal to 1, minimum is 0.
    SubSim = 1/n_sub * sum(SV.^2);
end