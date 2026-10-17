function [d_G] = GSPTB_ComputeOperatorDistance(G1,G2)

    % Distance between both operators
    d_G = norm(G1-G2,'fro')/norm(G2,'fro');
end