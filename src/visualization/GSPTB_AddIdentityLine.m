function [h] = GSPTB_AddIdentityLine(ax,Color,LineWidth,LineStyle)

    % Validates axes
    if isempty(ax) || ~isgraphics(ax,'axes')
        error('GSPTB:InvalidAxes', ...
            'ax must be a valid axes handle!');
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

    % Validates line width
    if ~isnumeric(LineWidth) || ~isscalar(LineWidth) || ...
            ~isfinite(LineWidth) || LineWidth <= 0
        error('GSPTB:InvalidLineWidth', ...
            'LineWidth must be a positive finite scalar!');
    end

    % Validates line style
    LineStyle = string(LineStyle);
    ValidStyles = ["-","--",":","-."];

    if ~isscalar(LineStyle) || ~ismember(LineStyle,ValidStyles)
        error('GSPTB:InvalidLineStyle', ...
            'LineStyle must be ''-'', ''--'', '':'', or ''-.''!');
    end

    % Determines the common visible range
    XLimits = xlim(ax);
    YLimits = ylim(ax);

    LowerLimit = max(XLimits(1),YLimits(1));
    UpperLimit = min(XLimits(2),YLimits(2));

    if LowerLimit >= UpperLimit
        error('GSPTB:IncompatibleAxisLimits', ...
            'The x- and y-axis limits do not overlap!');
    end

    % Preserve hold state
    WasHeld = ishold(ax);
    hold(ax,'on');

    % Add identity line
    hLine = line(ax, ...
        [LowerLimit UpperLimit], ...
        [LowerLimit UpperLimit], ...
        'Color',Color, ...
        'LineWidth',LineWidth, ...
        'LineStyle',LineStyle);

    % Restore hold state
    if ~WasHeld
        hold(ax,'off');
    end

    % Return handles
    h.Axes = ax;
    h.Line = hLine;
end