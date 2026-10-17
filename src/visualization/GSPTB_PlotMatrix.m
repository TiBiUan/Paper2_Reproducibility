% ax is the handle of the axes on which to plot
% M is a 2D array with the data to plot
% Colormap is an N x 3 color map
% ColorLimits specifies the color range; set to min max by default
function [h] = GSPTB_PlotMatrix(ax,M,ColorMap,ColorLimits)

    % If we provide an empty axis handle, then we create what we need
    if isempty(ax)
        f = figure;
        ax = axes(f);

    % Else, we only get the underlying figure handle from what was given as
    % input axis
    else
        f = ancestor(ax,'figure');
    end

    % Verifications for inputs

    % Empty?
    if isempty(M)
        error('GSPTB:EmptyInput','M must be non-empty!');
    end
    
    % Proper matrix given?
    if ~isnumeric(M) || ~ismatrix(M)
        error('GSPTB:InvalidMatrix', ...
            'M must be a numeric two-dimensional matrix!');
    end
    
    % Proper colormap fed?
    if ~isnumeric(ColorMap) || isempty(ColorMap) || ~ismatrix(ColorMap) || size(ColorMap,2) ~= 3
        error('GSPTB:InvalidColorMap', ...
            'ColorMap must be a non-empty numeric N-by-3 array!');
    end

    if any(~isfinite(ColorMap(:))) || any(ColorMap(:) < 0 | ColorMap(:) > 1)
        error('GSPTB:InvalidColorMap', ...
            'ColorMap values must be finite and lie between 0 and 1!');
    end
    
    % Proper limits for colors given?
    if ~isempty(ColorLimits) && (~isnumeric(ColorLimits) || ...
            numel(ColorLimits) ~= 2 || any(~isfinite(ColorLimits)) || ...
         ColorLimits(1) >= ColorLimits(2))
        error('GSPTB:InvalidColorLimits', ...
            'ColorLimits must contain two increasing numeric values!');
    end

    % Adjusts inputs
    ColorLimits = ColorLimits(:);

    % Preserve the initial hold state
    WasHeld = ishold(ax);
    hold(ax,'on');

    % Plotting per se
    hImage = imagesc(ax,M);

    % Restore the original hold state
    if ~WasHeld
        hold(ax,'off');
    end

    % If we have NaN entries, sets as transparent
    hImage.AlphaData = ~isnan(M);

    axis(ax,'tight');
    set(ax,'YDir','Reverse','Box','off');
    
    % Coloring
    colormap(ax,ColorMap);
    
    if ~isempty(ColorLimits)
        clim(ax,ColorLimits);
    end

    % Preparing output
    h.Figure = f;
    h.Axes = ax;
    h.Image = hImage;
end