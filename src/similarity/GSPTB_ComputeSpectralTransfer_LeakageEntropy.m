function [LeakageEntropy] = GSPTB_ComputeSpectralTransfer_LeakageEntropy(T)
    
    % Number of harmonics
    N = size(T,1);

    % Initialization
    LeakageEntropy = NaN(N,1);

    % Computes leakage entropy
    for j = 1:N

        % Samples the column of T denoting how harmonic j from U2 relates
        % to the basis U1
        t = T(:,j);
        t = t(t > 0);

        % Leakage entropy gives insight into the resulting distribution
        LeakageEntropy(j) = -sum(t.*log(t));
    end
end