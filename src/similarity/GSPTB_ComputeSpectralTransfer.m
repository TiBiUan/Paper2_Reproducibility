function [T] = GSPTB_ComputeSpectralTransfer(U1,U2)

    if ~isequal(size(U1),size(U2))
        error('U1 and U2 must have the same dimensions!');
    end
    
    T = abs(U1'*U2).^2;
end