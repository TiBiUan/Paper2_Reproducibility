% ax is the handle of the axes on which to plot
% Time has size 1 x n_timepoints
% Data has size n_entities x n_timepoints
% Colormap is an N x 3 color map
% ColorLimits specifies the color range; set to min max by default
function h = GSPTB_PlotTemporalMatrix(ax,Time,Data,ColorMap,ColorLimits)

    % Create axes if needed
    if isempty(ax)
        f = figure;
        ax = axes(f);
    else
        f = ancestor(ax,'figure');
    end

    % Validate inputs
    if isempty(Time) || isempty(Data)
        error('GSPTB:EmptyInput', ...
            'Time and Data must be non-empty!');
    end

    if ~isvector(Time)
        error('GSPTB:InvalidTime', ...
            'Time must be a vector!');
    end

    if any(~isfinite(Time))
        error('GSPTB:InvalidTime', ...
            'Time values must be finite!');
    end

    if any(diff(Time) <= 0)
        error('GSPTB:InvalidTime', ...
            'Time values must be strictly increasing!');
    end

    if ~isnumeric(Data) || ~ismatrix(Data)
        error('GSPTB:InvalidData', ...
            'Data must be a numeric two-dimensional matrix!');
    end

    % Plotting does not handle uneven spacing between time points, so we
    % explicitly check whether the provided time coordinates are equally
    % spaced overall
    dTime = diff(Time);

    if numel(dTime) > 1 && max(abs(dTime-mean(dTime))) > 1e-10*max(1,abs(mean(dTime)))
        warning('GSPTB:IrregularTimeSampling', ...
            ['Time is not regularly sampled. imagesc renders ' ...
             'columns with equal widths!']);
    end

    % Processing inputs
    Time = Time(:);

    if size(Data,2) ~= numel(Time)
        error('GSPTB:DimensionMismatch', ...
            ['Data has %d columns, but Time contains ' ...
             '%d elements!'], ...
            size(Data,2),numel(Time));
    end

    % Preserve the initial hold state
    WasHeld = ishold(ax);
    hold(ax,'on');

    % Plotting
    hImage = imagesc(ax,Time,1:size(Data,1),Data);

    % Restore the original hold state
    if ~WasHeld
        hold(ax,'off');
    end

    axis(ax,'tight');
    set(ax,'YDir','reverse','Box','off');

    colormap(ax,ColorMap);

    if ~isempty(ColorLimits)
        clim(ax,ColorLimits);
    end

    % Return handles
    h.Figure = f;
    h.Axes = ax;
    h.Image = hImage;
end