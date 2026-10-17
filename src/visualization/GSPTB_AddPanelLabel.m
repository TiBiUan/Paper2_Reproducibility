function [h] = GSPTB_AddPanelLabel(ax,Label,Position,FontSize)

    % Validates axes
    if isempty(ax) || ~isgraphics(ax,'axes')
        error('GSPTB:InvalidAxes', ...
            'ax must be a valid axes handle!');
    end

    % Validates label
    if ischar(Label)
        Label = string(Label);
    end

    if ~isstring(Label) || ~isscalar(Label) || ...
            ismissing(Label) || strlength(Label) == 0
        error('GSPTB:InvalidLabel', ...
            'Label must be a non-empty text scalar!');
    end

    % Validates normalized position
    if ~isnumeric(Position) || ~isreal(Position) || ...
            ~isequal(size(Position),[1 2]) || ...
            any(~isfinite(Position))
        error('GSPTB:InvalidPosition', ...
            'Position must be a finite 1-by-2 numeric vector!');
    end

    % Validates font size
    if ~isnumeric(FontSize) || ~isscalar(FontSize) || ...
            ~isfinite(FontSize) || FontSize <= 0
        error('GSPTB:InvalidFontSize', ...
            'FontSize must be a positive finite scalar!');
    end

    % Adds label in normalized axes coordinates
    hText = text(ax,Position(1),Position(2),Label, ...
        'Units','normalized', ...
        'FontSize',FontSize, ...
        'FontWeight','bold', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'Clipping','off');

    % Returns handles
    h.Axes = ax;
    h.Text = hText;
end