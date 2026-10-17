function [Err] = GSPTB_ComputeProjectionError(T,Windows)

    % Number of harmonics at hand
    N = size(T,1);
    
    % Checks on inputs
    if size(T,2) ~= N
        error('T must be a square transfer matrix!');
    end
    
    if any(Windows < 0) || any(mod(Windows,1) ~= 0)
        error('Windows must contain only non-negative integer values!');
    end
    
    % Projection errors initialization
    Err = zeros(N,numel(Windows));
    
    % We examine all the considered window schemes
    for nw = 1:numel(Windows)
    
        % We sample the current window "half-size" (this is "how many
        % spectral contributions we sample on one side of where we lie)
        w = Windows(nw);
    
        % For each harmonic from U2, we compute the error
        for j = 1:N

            % These are the indices from U1 to consider, centered at
            % harmonic j from U2
            idx = max(1,j-w):min(N,j+w);

            % This is the error; thanks to orthonormality of the bases, we
            % could simplify the expression to depend on T
            Err(j,nw) = 1 - sum(T(idx,j));
        end
    end
end