function [h] = GSPTB_AddVerticalReference(ax,X,Color,LineWidth,LineStyle)

    % Validates axes
    if isempty(ax) || ~isgraphics(ax,'axes')
        error('GSPTB:InvalidAxes', ...
            'ax must be a valid axes handle!');
    end

    % Validates reference positions
    if isempty(X) || ~isnumeric(X) || ~isvector(X) || ...
            any(~isfinite(X))
        error('GSPTB:InvalidReferencePosition', ...
            'X must be a non-empty finite numeric vector!');
    end

    X = X(:);

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
    ValidStyles = ["-","--",":","-."];

    LineStyle = string(LineStyle);

    if ~isscalar(LineStyle) || ...
            ~ismember(LineStyle,ValidStyles)
        error('GSPTB:InvalidLineStyle', ...
            'LineStyle must be ''-'', ''--'', '':'', or ''-.''!');
    end

    % Preserve hold state
    WasHeld = ishold(ax);
    hold(ax,'on');

    % Plot references
    N = numel(X);
    hLines = gobjects(N,1);

    for n = 1:N
        hLines(n) = xline(ax,X(n), ...
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