function [RequiredSpectralWindow] = ...
    GSPTB_ComputeProjectionError_RequiredSpectralWindow(Err,Windows,Tolerance)

    if size(Err,2) ~= numel(Windows)
        error('Number of columns in Err must match the number of window sizes!');
    end
    
    if any(diff(Windows) <= 0)
        error('Windows must contain strictly increasing values!');
    end

    % Err has size N (number of assessed harmonics) x W (number of windows)
    N = size(Err,1);

    % Initializes
    RequiredSpectralWindow = NaN(N,1);
    
    for j = 1:N

        % Finds the first index for which the error becomes lower than the
        % tolerance
        id = find(Err(j,:) <= Tolerance,1,'first');

        if ~isempty(id)
            RequiredSpectralWindow(j) = Windows(id);
        end
    end
end