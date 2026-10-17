function [MU,STD,VAR] = GSPTB_ComputeStatisticsFromSS(S,SS,N)

    MU = S./N;
    VAR = (SS-S.^2./N)./(N-1);

    if isscalar(N)
        if N == 0
            MU(:) = NaN;
            VAR(:) = NaN;
        elseif N == 1
            VAR(:) = NaN;
        else
            VAR = max(VAR,0);
        end
    else
        MU(N == 0) = NaN;

        ValidVar = N >= 2;
        VAR(ValidVar) = max(VAR(ValidVar),0);
        VAR(~ValidVar) = NaN;
    end

    STD = sqrt(VAR);
end