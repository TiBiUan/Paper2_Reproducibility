function [h] = GSPTB_PlotBoxDistribution( ...
    ax,Values,Groups,BoxColors,varargin)

    % GSPTB_PlotBoxDistribution
    %
    % Displays grouped numeric observations using box charts, with
    % optional jittered raw-data points and customizable x-axis labels.
    %
    % INPUTS
    %   ax
    %       Destination axes. Pass [] to create a new figure and axes.
    %
    %   Values
    %       Numeric vector containing the observations.
    %
    %   Groups
    %       Numeric, categorical, or string vector assigning each
    %       observation to a group.
    %
    %   BoxColors
    %       Either one RGB triplet or one RGB triplet per group.
    %
    % NAME-VALUE OPTIONS
    %   'ShowPoints'
    %       Logical scalar indicating whether individual observations
    %       should be displayed. Default: false.
    %
    %   'PointSize'
    %       Positive scalar controlling marker size. Default: 10.
    %
    %   'PointAlpha'
    %       Scalar between 0 and 1 controlling point transparency.
    %       Default: 0.20.
    %
    %   'JitterWidth'
    %       Non-negative horizontal jitter half-width. Default: 0.12.
    %
    %   'MaximumPointsPerGroup'
    %       Maximum number of raw observations displayed in each group.
    %       The box chart always uses all observations. Use Inf to display
    %       every point. Default: Inf.
    %
    %   'JitterSeed'
    %       Non-negative integer used for reproducible point placement and
    %       optional display subsampling. Default: 1.
    %
    %   'XTickPositions'
    %       Numeric vector specifying custom x-axis tick positions.
    %       Default: 1:G.
    %
    %   'XTickLabels'
    %       Labels corresponding to XTickPositions. By default, the unique
    %       group values are shown.
    %
    % OUTPUT
    %   h
    %       Structure containing figure, axes, box-chart, and scatter
    %       handles, together with group and plotting metadata.


    %% Parse optional inputs

    Parser = inputParser;

    addParameter(Parser,'ShowPoints',false, ...
        @(x) islogical(x) && isscalar(x));

    addParameter(Parser,'PointSize',10, ...
        @(x) isnumeric(x) && isreal(x) && isscalar(x) && ...
             isfinite(x) && x > 0);

    addParameter(Parser,'PointAlpha',0.20, ...
        @(x) isnumeric(x) && isreal(x) && isscalar(x) && ...
             isfinite(x) && x >= 0 && x <= 1);

    addParameter(Parser,'JitterWidth',0.12, ...
        @(x) isnumeric(x) && isreal(x) && isscalar(x) && ...
             isfinite(x) && x >= 0);

    addParameter(Parser,'MaximumPointsPerGroup',Inf, ...
        @(x) isnumeric(x) && isreal(x) && isscalar(x) && ...
             x > 0 && (isinf(x) || ...
             (isfinite(x) && mod(x,1) == 0)));

    addParameter(Parser,'JitterSeed',1, ...
        @(x) isnumeric(x) && isreal(x) && isscalar(x) && ...
             isfinite(x) && x >= 0 && mod(x,1) == 0);

    addParameter(Parser,'XTickPositions',[], ...
        @(x) isempty(x) || ...
             (isnumeric(x) && isreal(x) && isvector(x) && ...
              all(isfinite(x))));

    addParameter(Parser,'XTickLabels',[], ...
        @(x) isempty(x) || isstring(x) || iscellstr(x) || ...
             iscategorical(x) || isnumeric(x));

    parse(Parser,varargin{:});
    Options = Parser.Results;


    %% Create or validate axes

    if isempty(ax)

        f = figure;
        ax = axes(f);

    elseif ~isscalar(ax) || ~isgraphics(ax,'axes')

        error('GSPTB:InvalidAxes', ...
            'ax must be empty or a valid scalar axes handle!');

    else

        f = ancestor(ax,'figure');

    end


    %% Validate Values and Groups

    if isempty(Values) || isempty(Groups)

        error('GSPTB:EmptyInput', ...
            'Groups and Values must be non-empty!');

    end

    if ~isnumeric(Values) || ~isreal(Values) || ~isvector(Values)

        error('GSPTB:InvalidData', ...
            'Values must be a real numeric vector!');

    end

    if ~(isnumeric(Groups) || iscategorical(Groups) || ...
            isstring(Groups)) || ~isvector(Groups)

        error('GSPTB:InvalidGroups', ...
            ['Groups must be a numeric, categorical, ' ...
             'or string vector!']);

    end


    %% Vectorize inputs

    Values = Values(:);
    Groups = Groups(:);

    if numel(Values) ~= numel(Groups)

        error('GSPTB:DimensionMismatch', ...
            ['Values contains %d elements, but Groups contains ' ...
             '%d elements!'],numel(Values),numel(Groups));

    end


    %% Remove paired missing observations

    if isnumeric(Groups)

        Valid = ~isnan(Values) & ~isnan(Groups);

    else

        Valid = ~isnan(Values) & ~ismissing(Groups);

    end

    Values = Values(Valid);
    Groups = Groups(Valid);

    if isempty(Values)

        error('GSPTB:EmptyInput', ...
            ['No valid observations remain after removing ' ...
             'missing data!']);

    end

    if any(~isfinite(Values))

        error('GSPTB:InvalidData', ...
            'Values must contain only finite observations!');

    end


    %% Determine group ordering

    UniqueGroups = unique(Groups,'stable');
    G = numel(UniqueGroups);


    %% Validate box colors

    if ~isnumeric(BoxColors) || ~isreal(BoxColors) || ...
            ~ismatrix(BoxColors) || size(BoxColors,2) ~= 3 || ...
            ~(size(BoxColors,1) == 1 || size(BoxColors,1) == G) || ...
            any(~isfinite(BoxColors(:))) || ...
            any(BoxColors(:) < 0 | BoxColors(:) > 1)

        error('GSPTB:InvalidColor', ...
            ['BoxColors must be a valid 1-by-3 RGB triplet or ' ...
             'a G-by-3 array containing one color per group!']);

    end

    if size(BoxColors,1) == 1
        BoxColors = repmat(BoxColors,G,1);
    end


    %% Configure custom tick positions and labels

    if isempty(Options.XTickPositions)
        TickPositions = 1:G;
    else
        TickPositions = Options.XTickPositions(:).';
    end

    if isempty(Options.XTickLabels)

        if numel(TickPositions) ~= G

            error('GSPTB:TickLabelMismatch', ...
                ['Custom XTickPositions require corresponding ' ...
                 'XTickLabels unless their number equals the ' ...
                 'number of groups.']);

        end

        TickLabels = string(UniqueGroups);

    else

        TickLabels = string(Options.XTickLabels(:));

        if numel(TickLabels) ~= numel(TickPositions)

            error('GSPTB:TickLabelMismatch', ...
                ['XTickLabels contains %d labels, whereas ' ...
                 'XTickPositions contains %d positions!'], ...
                 numel(TickLabels),numel(TickPositions));

        end

    end


    %% Preserve the original hold state

    WasHeld = ishold(ax);
    hold(ax,'on');


    %% Prepare output handles

    hBox = gobjects(G,1);
    hPoints = gobjects(G,1);

    DisplayedPointIndices = cell(G,1);


    %% Preserve the caller's random-number state

    PreviousRandomState = rng;
    RestoreRandomState = onCleanup(@() rng(PreviousRandomState));

    rng(Options.JitterSeed,'twister');


    %% Plot each distribution

    for g = 1:G

        Mask = Groups == UniqueGroups(g);
        CurrentValues = Values(Mask);
        NumberOfValues = numel(CurrentValues);

        % Optionally display raw observations.
        if Options.ShowPoints

            NumberToDisplay = min( ...
                NumberOfValues,Options.MaximumPointsPerGroup);

            if NumberToDisplay < NumberOfValues

                SelectedIndices = randperm( ...
                    NumberOfValues,NumberToDisplay).';

            else

                SelectedIndices = (1:NumberOfValues).';

            end

            DisplayedValues = CurrentValues(SelectedIndices);

            XCoordinates = repmat(g,NumberToDisplay,1);

            if Options.JitterWidth > 0

                XCoordinates = XCoordinates + ...
                    (2*rand(NumberToDisplay,1)-1) * ...
                    Options.JitterWidth;

            end

            hPoints(g) = scatter(ax, ...
                XCoordinates, ...
                DisplayedValues, ...
                Options.PointSize, ...
                BoxColors(g,:), ...
                'filled', ...
                'MarkerFaceAlpha',Options.PointAlpha, ...
                'MarkerEdgeAlpha',Options.PointAlpha);

            DisplayedPointIndices{g} = SelectedIndices;

        else

            DisplayedPointIndices{g} = zeros(0,1);

        end

        % The box chart always uses every valid observation.
        hBox(g) = boxchart(ax, ...
            repmat(g,NumberOfValues,1), ...
            CurrentValues, ...
            'BoxFaceColor',BoxColors(g,:), ...
            'MarkerStyle','none');

    end


    %% Restore original hold state

    if ~WasHeld
        hold(ax,'off');
    end


    %% Apply axes aesthetics

    xticks(ax,TickPositions);
    xticklabels(ax,TickLabels);

    set(ax,'Box','off');


    %% Return handles and metadata

    h.Figure = f;
    h.Axes = ax;
    h.Boxes = hBox;
    h.Points = hPoints;
    h.Values = Values;
    h.Groups = Groups;
    h.UniqueGroups = UniqueGroups;
    h.BoxColors = BoxColors;
    h.TickPositions = TickPositions;
    h.TickLabels = TickLabels;
    h.DisplayedPointIndices = DisplayedPointIndices;
    h.Options = Options;
end