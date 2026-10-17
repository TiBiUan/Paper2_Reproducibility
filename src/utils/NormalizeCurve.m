function [y2] = NormalizeCurve(y)

    y2 = y - min(y);
    y2 = y2/max(y2);
end