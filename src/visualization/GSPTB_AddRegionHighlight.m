function [h] = GSPTB_AddRegionHighlight(ax,Range,Orientation,Color,Alpha)

    % Validates axes
    if isempty(ax) || ~isgraphics(ax,'axes')
        error('GSPTB:InvalidAxes', ...
            'ax must be a valid axes handle!');
    end

    % Validates range
    if ~isnumeric(Range) || ...
            ~isequal(size(Range),[1 2]) || ...
            any(~isfinite(Range)) || ...
            Range(1) >= Range(2)
        error('GSPTB:InvalidRange', ...
            ['Range must contain two finite values ' ...
             'in strictly increasing order!']);
    end

    % Validates orientation
    Orientation = lower(string(Orientation));

    if ~isscalar(Orientation) || ...
            ~ismember(Orientation,["vertical","horizontal"])
        error('GSPTB:InvalidOrientation', ...
            'Orientation must be ''vertical'' or ''horizontal''!');
    end

    % Validates color
    if ~isnumeric(Color) || ...
            ~isequal(size(Color),[1 3]) || ...
            any(~isfinite(Color)) || ...
            any(Color < 0 | Color > 1)
        error('GSPTB:InvalidColor', ...
            ['Color must be a finite 1-by-3 RGB triplet ' ...
             'with values between 0 and 1!']);
    end

    % Validates transparency
    if ~isnumeric(Alpha) || ...
            ~isscalar(Alpha) || ...
            ~isfinite(Alpha) || ...
            Alpha < 0 || Alpha > 1
        error('GSPTB:InvalidAlpha', ...
            'Alpha must be a scalar between 0 and 1!');
    end

    % Preserves hold state
    WasHeld = ishold(ax);
    hold(ax,'on');

    % Current axis limits
    XLimits = xlim(ax);
    YLimits = ylim(ax);

    % Draws patch
    if Orientation == "vertical"
        X = [Range(1) Range(2) Range(2) Range(1)];
        Y = [YLimits(1) YLimits(1) YLimits(2) YLimits(2)];
    else
        X = [XLimits(1) XLimits(2) XLimits(2) XLimits(1)];
        Y = [Range(1) Range(1) Range(2) Range(2)];
    end

    % Plots the patch
    hPatch = patch(ax,X,Y,Color,...
        'EdgeColor','none',...
        'FaceAlpha',Alpha);

    % Keeps plotted data visible
    uistack(hPatch,'bottom');

    % Restores hold state
    if ~WasHeld
        hold(ax,'off');
    end

    % Returns handles
    h.Axes = ax;
    h.Patch = hPatch;
end