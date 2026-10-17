function [A_Pop] = GSPTB_ComputePopulationAverageNonZero(A,idx)

    if isempty(idx)
        error('Vector idx should not be empty!');
    end

    A_Selected = cat(3,A{idx});

    % How many weights are non-zero
    C = sum(A_Selected > 0,3);

    % We average over non-null weights only, and make sure that whatever
    % entry is null in all subjects is zero in the output
    A_Pop = sum(A_Selected,3)./C;
    A_Pop(C == 0) = 0;
end