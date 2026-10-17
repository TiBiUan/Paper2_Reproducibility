function [D] = GSPTB_ComputeDifferenceMatrix(Lambda1,Lambda2)

    % Number of harmonics
    N = length(Lambda1);
        
    D = NaN(N,N);

    for i = 1:N
        for j = 1:N
            D(i,j) = abs(Lambda1(i)-Lambda2(j));
        end
    end
end