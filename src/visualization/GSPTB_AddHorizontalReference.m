function [h] = GSPTB_AddHorizontalReference( ...
    ax,Y,Color,LineWidth,LineStyle)

    % Validate axes
    if isempty(ax) || ~isgraphics(ax,'axes')
        error('GSPTB:InvalidAxes', ...
            'ax must be a valid axes handle!');
    end

    % Validate positions
    if isempty(Y) || ~isnumeric(Y) || ~isvector(Y) || ...
            any(~isfinite(Y))
        error('GSPTB:InvalidReferencePosition', ...
            'Y must be a non-empty finite numeric vector!');
    end

    Y = Y(:);

    % Validate color
    if ~isnumeric(Color) || ...
            ~isequal(size(Color),[1 3]) || ...
            any(~isfinite(Color)) || ...
            any(Color < 0 | Color > 1)
        error('GSPTB:InvalidColor', ...
            ['Color must be a finite 1-by-3 RGB triplet ' ...
             'with values between 0 and 1!']);
    end

    % Validate line width
    if ~isnumeric(LineWidth) || ~isscalar(LineWidth) || ...
            ~isfinite(LineWidth) || LineWidth <= 0
        error('GSPTB:InvalidLineWidth', ...
            'LineWidth must be a positive finite scalar!');
    end

    % Validate line style
    LineStyle = string(LineStyle);
    ValidStyles = ["-","--",":","-."];

    if ~isscalar(LineStyle) || ~ismember(LineStyle,ValidStyles)
        error('GSPTB:InvalidLineStyle', ...
            'LineStyle must be ''-'', ''--'', '':'', or ''-.''!');
    end

    % Preserve hold state
    WasHeld = ishold(ax);
    hold(ax,'on');

    % Plot references
    N = numel(Y);
    hLines = gobjects(N,1);

    for n = 1:N
        hLines(n) = yline(ax,Y(n), ...
            'Color',Color, ...
            'LineWidth',LineWidth, ...
            'LineStyle',LineStyle);

        hLines(n).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end

    % Restore hold state
    if ~WasHeld
        hold(ax,'off');
    end

    

    % Return handles
    h.Axes = ax;
    h.Lines = hLines;
end