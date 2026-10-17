function [A_Pop] = GSPTB_ComputePopulationAverage(A,idx)

    if isempty(idx)
        error('Vector idx should not be empty!');
    end

    A_Selected = cat(3,A{idx});

    % Very concise approach to averaging
    A_Pop = mean(A_Selected,3);
end