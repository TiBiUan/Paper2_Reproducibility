function [h] = GSPTB_PlotRelationship(ax,X,Y,PointColors,PointSize)
    
    % Create axes if needed
    if isempty(ax)
        f = figure;
        ax = axes(f);
    else
        f = ancestor(ax,'figure');
    end
    
    % Handling of data-related issues
    if isempty(X) || isempty(Y)
        error('GSPTB:EmptyInput', ...
            'X and Y must be non-empty!');
    end

    if ~isnumeric(X) || ~isnumeric(Y)
        error('GSPTB:InvalidData', ...
            'X and Y must be numeric arrays!');
    end

    if any(isinf(X)) || any(isinf(Y))
        error('GSPTB:InvalidData', ...
            'X and Y must not contain infinite values!');
    end

    X = X(:);
    Y = Y(:);

    if numel(X) ~= numel(Y)
        error('GSPTB:DimensionMismatch', ...
            ['X has %d elements, but Y contains ' ...
             '%d elements!'], ...
            numel(X),numel(Y));
    end

    % Handling of color input
    N = numel(X);

    if ~isnumeric(PointColors)
        error('GSPTB:InvalidPointColor', ...
            'PointColors must be numeric!');
    end

    % Checks regarding point colors, which must be of size 1 x 3 or N x 1
    if isequal(size(PointColors),[1 3])
    
        if any(~isfinite(PointColors)) || any(PointColors < 0 | PointColors > 1)
            error('GSPTB:InvalidPointColors', ...
                'RGB values must lie between 0 and 1!');
        end

    elseif iscolumn(PointColors) && numel(PointColors) == N
    
        if any(~isfinite(PointColors))
            error('GSPTB:InvalidPointColors','PointColors contains non-finite values!');
        end
    
    else
        error('GSPTB:InvalidPointColors', ...
            ['PointColors must be either a 1-by-3 RGB triplet ' ...
             'or an N-by-1 vector containing one scalar value per point!']);
    end

    % Checks regarding point sizes
    if ~isnumeric(PointSize) || isempty(PointSize) || ...
        any(~isfinite(PointSize(:))) || any(PointSize(:) <= 0)
        error('GSPTB:InvalidPointSize','PointSize must contain positive finite numeric values!');
    end
    
    if ~isscalar(PointSize) && numel(PointSize) ~= N
        error('GSPTB:DimensionMismatch',['PointSize must be scalar or contain one value ' ...
             'for each of the %d points!'],N);
    end
    
    PointSize = PointSize(:).';
           
    % Preserve the initial hold state
    WasHeld = ishold(ax);
    hold(ax,'on');

    % Plotting per se, always done the same
    hScatter = scatter(ax,X,Y,PointSize,PointColors,'filled');
    
    % Restore the original hold state
    if ~WasHeld
        hold(ax,'off');
    end

    % Further aesthetics
    set(ax,'Box','off');

    % Output to return
    h.Figure = f;
    h.Axes = ax;
    h.Scatter = hScatter;
end