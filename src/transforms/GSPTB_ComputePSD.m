function [PSD] = GSPTB_ComputePSD(X_hat)

    % Power spectral density, averaged across time points
    PSD = mean(abs(X_hat).^2,2);
end