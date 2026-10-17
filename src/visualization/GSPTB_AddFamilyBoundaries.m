function [h] = GSPTB_AddFamilyBoundaries(ax,Boundaries,Orientation,Color,LineWidth,LineStyle)

    % Validation for Axes
    if isempty(ax) || ~isgraphics(ax,'axes')
        error('GSPTB:InvalidAxes', ...
            'ax must be a valid axes handle!');
    end

    % Validation for Boundaries
    if isempty(Boundaries) || ~isnumeric(Boundaries) || ...
            ~isvector(Boundaries) || any(~isfinite(Boundaries))
        error('GSPTB:InvalidBoundaries', ...
            'Boundaries must be a non-empty finite numeric vector!');
    end

    Boundaries = Boundaries(:);

    % Validation for Orientation
    Orientation = lower(string(Orientation));

    if ~isscalar(Orientation) || ...
            ~ismember(Orientation,["vertical","horizontal","both"])
        error('GSPTB:InvalidOrientation', ...
            ['Orientation must be ''vertical'', ' ...
             '''horizontal'', or ''both''!']);
    end

    h.Axes = ax;
    h.Vertical = [];
    h.Horizontal = [];

    % We make use of our simpler helper functions
    if Orientation == "vertical" || Orientation == "both"
        h.Vertical = GSPTB_AddVerticalReference(ax,Boundaries,Color,LineWidth,LineStyle);
    end

    if Orientation == "horizontal" || Orientation == "both"
        h.Horizontal = GSPTB_AddHorizontalReference(ax,Boundaries,Color,LineWidth,LineStyle);
    end
end