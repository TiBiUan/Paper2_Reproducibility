function [SubspaceSimilarity,Idx_SubspaceSimilarity] = ...
    GSPTB_ComputeSubspaceSimilarityAcrossWindows(U1,U2,Windows,SimThresh)

    if ~isequal(size(U1),size(U2))
        error('U1 and U2 must have the same dimensions!');
    end
    
    if any(Windows < 0) || any(mod(Windows,1) ~= 0)
        error('Windows must contain only non-negative integer values!');
    end
    
    if any(diff(Windows) <= 0)
        error('Windows must contain strictly increasing values!');
    end
    
    if SimThresh < 0 || SimThresh > 1
        error('SimThresh must lie between 0 and 1!');
    end

    N = size(U1,2);

    SubspaceSimilarity = NaN(N,numel(Windows));
    Idx_SubspaceSimilarity = NaN(N,1);

    % Squared full overlap matrix
    C2 = (U1' * U2).^2;

    % 2-D cumulative sum, padded with zeros
    CS = zeros(N+1,N+1);
    CS(2:end,2:end) = cumsum(cumsum(C2,1),2);

    for j = 1:N

        for nw = 1:numel(Windows)

            w = Windows(nw);

            a = max(1,j-w);
            b = min(N,j+w);

            n_sub = b-a+1;

            % Sum of C2(a:b,a:b), obtained from the 2-D cumulative sum
            BlockSum = ...
                CS(b+1,b+1) ...
                - CS(a,b+1) ...
                - CS(b+1,a) ...
                + CS(a,a);

            SubspaceSimilarity(j,nw) = BlockSum / n_sub;
        end

        id = find(SubspaceSimilarity(j,:) >= SimThresh,1,'first');

        if ~isempty(id)
            Idx_SubspaceSimilarity(j) = Windows(id);
        end
    end
end