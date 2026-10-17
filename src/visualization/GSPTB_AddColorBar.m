function [h] = GSPTB_AddColorBar(ax,Label)

    % Validates axes
    if isempty(ax) || ~isgraphics(ax,'axes')
        error('GSPTB:InvalidAxes', ...
            'ax must be a valid axes handle!');
    end

    % Validates label
    if ischar(Label)
        Label = string(Label);
    end

    if ~isstring(Label) || ~isscalar(Label) || ismissing(Label)
        error('GSPTB:InvalidLabel', ...
            'Label must be a valid text scalar!');
    end

    % Create colorbar
    hColorBar = colorbar(ax);

    % Add label when provided
    if strlength(Label) > 0
        hColorBar.Label.String = Label;
    end

    % Return handles
    h.Axes = ax;
    h.ColorBar = hColorBar;
end