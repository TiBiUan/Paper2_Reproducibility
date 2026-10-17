function [h] = GSPTB_PlotCategoricalMatrix(ax,Data,CategoryValues,CategoryColors)

    % Create axes if needed
    if isempty(ax)
        f = figure;
        ax = axes(f);
    else
        f = ancestor(ax,'figure');
    end

    % Validate Data
    if isempty(Data)
        error('GSPTB:EmptyInput', ...
            'Data must be non-empty!');
    end

    if ~isnumeric(Data) || ~ismatrix(Data)
        error('GSPTB:InvalidData', ...
            'Data must be a numeric two-dimensional matrix!');
    end

    if any(isinf(Data(:)))
        error('GSPTB:InvalidData', ...
            'Data must not contain infinite values!');
    end

    % Validate category values
    if ~isnumeric(CategoryValues) || ...
            ~isvector(CategoryValues) || ...
            isempty(CategoryValues) || ...
            any(~isfinite(CategoryValues))
        error('GSPTB:InvalidCategories', ...
            'CategoryValues must be a non-empty finite numeric vector!');
    end

    CategoryValues = CategoryValues(:);

    if numel(unique(CategoryValues)) ~= numel(CategoryValues)
        error('GSPTB:InvalidCategories', ...
            'CategoryValues must not contain repeated values!');
    end

    % Validate category colors
    K = numel(CategoryValues);

    if ~isnumeric(CategoryColors) || ...
            ~ismatrix(CategoryColors) || ...
            ~isequal(size(CategoryColors),[K 3]) || ...
            any(~isfinite(CategoryColors(:))) || ...
            any(CategoryColors(:) < 0 | CategoryColors(:) > 1)
        error('GSPTB:InvalidColor', ...
            ['CategoryColors must be a K-by-3 RGB array with ' ...
             'values between 0 and 1!']);
    end

    % Convert category values to colormap indices
    [IsValid,ColorIndex] = ismember(Data,CategoryValues);

    if any(~IsValid(:) & ~isnan(Data(:)))
        error('GSPTB:InvalidCategory', ...
            'Data contains values not listed in CategoryValues!');
    end

    % Preserve the initial hold state
    WasHeld = ishold(ax);
    hold(ax,'on');

    % Plot category indices
    hImage = imagesc(ax,ColorIndex);

    % Restore the original hold state
    if ~WasHeld
        hold(ax,'off');
    end

    % Missing entries are transparent
    hImage.AlphaData = ~isnan(Data);

    colormap(ax,CategoryColors);
    clim(ax,[0.5 K + 0.5]);

    axis(ax,'tight');
    set(ax,'YDir','reverse','Box','off');

    % Return handles and category metadata
    h.Figure = f;
    h.Axes = ax;
    h.Image = hImage;
    h.CategoryValues = CategoryValues;
    h.CategoryColors = CategoryColors;
end