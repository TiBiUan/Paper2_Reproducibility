function [K] = GSPTB_IsSameFamily(Parents)

    for s1 = 1:size(Parents,1)
        for s2 = 1:size(Parents,1)

            if sum(ismember(Parents(s1,:),Parents(s2,:)) > 0)
                K(s1,s2) = 1;
            else
                K(s1,s2) = 0;
            end
        end
    end
end