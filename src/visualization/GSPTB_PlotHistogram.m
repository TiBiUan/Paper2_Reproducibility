function [h] = GSPTB_PlotHistogram(ax,Values,Nbins,Normalization,Color,Alpha)

    % Create axes if needed
    if isempty(ax)
        f = figure;
        ax = axes(f);
    else
        f = ancestor(ax,'figure');
    end
    
    % Tests for Values
    if isempty(Values)
        error('GSPTB:EmptyInput', ...
            'Values must be non-empty!');
    end
    
    if ~isnumeric(Values) || ~isvector(Values)
        error('GSPTB:InvalidData', ...
            'Values must be a numeric vector!');
    end
    
    if any(isinf(Values))
        error('GSPTB:InvalidData', ...
            'Values must not contain infinite values!');
    end
    
    Values = Values(:);

    Values = Values(~isnan(Values));

    if isempty(Values)
        error('GSPTB:EmptyInput','Values contain no finite observations!');
    end

    % Tests for Nbins
    if ~isnumeric(Nbins) || ~isscalar(Nbins) || ~isfinite(Nbins) || ...
            Nbins < 1 || mod(Nbins,1) ~= 0
        error('GSPTB:InvalidBinCount','Nbins must be a positive integer scalar!');
    end

    Normalization = lower(string(Normalization));
    ValidNormalizations = ["pdf","cdf","count","probability"];

    if ~ismember(Normalization,ValidNormalizations)
        error('GSPTB:InvalidNormalization', ...
            'Normalization must be ''pdf'', ''cdf'', ''count'', or ''probability''!');
    end

    % Tests for Color
    if ~isnumeric(Color) || ~isequal(size(Color),[1 3]) || ...
        any(~isfinite(Color)) || any(Color < 0 | Color > 1)
        error('GSPTB:InvalidColor',...
            'Color must be a finite 1-by-3 RGB triplet between 0 and 1!');
    end

    % Tests for transparency
    if ~isnumeric(Alpha) || ~isscalar(Alpha) || ...
        ~isfinite(Alpha) || Alpha < 0 || Alpha > 1
        error('GSPTB:InvalidAlpha', ...
            'Alpha must be a scalar between 0 and 1!');
    end

    % Preserve the initial hold state
    WasHeld = ishold(ax);
    hold(ax,'on');

    % Plotting per se
    hHistogram = histogram(ax,Values,Nbins,'Normalization',Normalization, ...
    'FaceColor',Color,'FaceAlpha',Alpha,'EdgeColor','none');

    % Restore the original hold state
    if ~WasHeld
        hold(ax,'off');
    end

    % Further aesthetics
    set(ax,'Box','off');

    h.Figure = f;
    h.Axes = ax;
    h.Histogram = hHistogram;
end

