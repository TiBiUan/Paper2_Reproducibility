function [G] = GSPTB_ComputeGraphFilterOperator(U,g)

    if size(U,2) ~= numel(g)
        error('Length of g must equal the number of harmonics in U!');
    end

    % We create the operator linked to the frequency profile g
    G = U*diag(g)*U';
end