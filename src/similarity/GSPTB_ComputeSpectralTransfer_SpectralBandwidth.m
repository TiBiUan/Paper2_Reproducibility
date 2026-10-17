function [SpectralBandwidth] = GSPTB_ComputeSpectralTransfer_SpectralBandwidth(T,Lambda1,Lambda2,Type)
    
    % Number of harmonics
    N = size(T,1);

    % We build the "difference matrix", denoting how "distant" harmonic i
    % from basis 1 is from harmonic j from basis 2 in D(i,j). We can
    % consider indices, or eigenvalues
    switch Type
        case 'index'
            D = GSPTB_ComputeDifferenceMatrix(1:N,1:N);
        case 'eigenvalue'
            D = GSPTB_ComputeDifferenceMatrix(Lambda1,Lambda2);
        otherwise
            error('Unrecognized type of quantity to contrast!');
    end

    % Spectral bandwidth
    SpectralBandwidth = sum(D.*T,1)';
end