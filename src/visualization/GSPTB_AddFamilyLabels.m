function [h] = GSPTB_AddFamilyLabels(ax,FamilySizes,Labels,Location,FontSize,Color)

    % Validates axes
    if isempty(ax) || ~isgraphics(ax,'axes')
        error('GSPTB:InvalidAxes', ...
            'ax must be a valid axes handle!');
    end

    % Validates family sizes
    if ~isnumeric(FamilySizes) || ~isreal(FamilySizes) || ...
            ~isvector(FamilySizes) || isempty(FamilySizes) || ...
            any(~isfinite(FamilySizes)) || ...
            any(FamilySizes <= 0) || ...
            any(mod(FamilySizes,1) ~= 0)
        error('GSPTB:InvalidFamilySizes', ...
            ['FamilySizes must be a non-empty vector of ' ...
             'positive integer values!']);
    end

    FamilySizes = FamilySizes(:);
    NFamilies = numel(FamilySizes);

    % Validates labels
    if ischar(Labels) || iscellstr(Labels) || iscategorical(Labels)
        Labels = string(Labels);
    end

    if ~isstring(Labels) || ~isvector(Labels) || isempty(Labels)
        error('GSPTB:InvalidLabels', ...
            'Labels must be a non-empty string vector!');
    end

    Labels = Labels(:);

    if any(ismissing(Labels)) || any(strlength(Labels) == 0)
        error('GSPTB:InvalidLabels', ...
            'Labels must not contain missing or empty entries!');
    end

    if numel(Labels) ~= NFamilies
        error('GSPTB:DimensionMismatch', ...
            ['FamilySizes defines %d families, but Labels ' ...
             'contains %d entries!'], ...
            NFamilies,numel(Labels));
    end

    % Validates location
    Location = lower(string(Location));

    if ~isscalar(Location) || ...
            ~ismember(Location,["top","bottom"])
        error('GSPTB:InvalidLocation', ...
            'Location must be ''top'' or ''bottom''!');
    end

    % Validates font size
    if ~isnumeric(FontSize) || ~isreal(FontSize) || ...
            ~isscalar(FontSize) || ...
            ~isfinite(FontSize) || FontSize <= 0
        error('GSPTB:InvalidFontSize', ...
            'FontSize must be a positive finite scalar!');
    end

    % Validates color
    if ~isnumeric(Color) || ~isreal(Color) || ...
            ~isequal(size(Color),[1 3]) || ...
            any(~isfinite(Color)) || ...
            any(Color < 0 | Color > 1)
        error('GSPTB:InvalidColor', ...
            ['Color must be a finite 1-by-3 RGB triplet ' ...
             'with values between 0 and 1!']);
    end

    % Compute family ranges and centers
    FamilyEnds = cumsum(FamilySizes);
    FamilyStarts = [1; FamilyEnds(1:end-1) + 1];
    Centers = (FamilyStarts + FamilyEnds) / 2;

    % Determine vertical label position
    YLimits = ylim(ax);
    Offset = 0.02 * diff(YLimits);

    if Location == "top"
        YPosition = YLimits(2) + Offset;
        VerticalAlignment = 'bottom';
    else
        YPosition = YLimits(1) - Offset;
        VerticalAlignment = 'top';
    end

    % Preserves hold state
    WasHeld = ishold(ax);
    hold(ax,'on');

    % Adds labels
    hText = gobjects(NFamilies,1);

    % Adds the text labels
    for k = 1:NFamilies
        hText(k) = text(ax, ...
            Centers(k), ...
            YPosition, ...
            Labels(k), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment',VerticalAlignment, ...
            'FontSize',FontSize, ...
            'Color',Color, ...
            'Clipping','off');
    end

    % Restores hold state
    if ~WasHeld
        hold(ax,'off');
    end

    % Returns handles and derived positions
    h.Axes = ax;
    h.Text = hText;
    h.Centers = Centers;
    h.FamilyStarts = FamilyStarts;
    h.FamilyEnds = FamilyEnds;
end