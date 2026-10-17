function [G] = GSPTB_CreateRectangularFilter(id,N)

    G = zeros(N,N);
    G(sub2ind([N,N],id,id)) = 1;
end