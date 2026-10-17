function [h] = GSPTB_PlotGraph(ax,Adjacency,Coordinates,NodeColors,...
    NodeSizes,EdgeColor,EdgeWidthRange)

    % Create axes if needed
    if isempty(ax)
        f = figure;
        ax = axes(f);
    else
        f = ancestor(ax,'figure');
    end

    % Validations on the Adjacency and Coordinates inputs
    if isempty(Adjacency) || isempty(Coordinates)
        error('GSPTB:EmptyInput', ...
            'Adjacency and Coordinates must be non-empty!');
    end
    
    if ~isnumeric(Adjacency) || ~isreal(Adjacency) || ...
        ~ismatrix(Adjacency) || ...
        size(Adjacency,1) ~= size(Adjacency,2)
        error('GSPTB:InvalidAdjacency', ...
            'Adjacency must be a square numeric matrix of real numbers!');
    end
    
    if any(~isfinite(Adjacency(:)))
        error('GSPTB:InvalidAdjacency', ...
            'Adjacency must contain only finite values!');
    end
    
    N = size(Adjacency,1);
    
    if ~isnumeric(Coordinates) || ~ismatrix(Coordinates) || ...
            size(Coordinates,1) ~= N || ~isreal(Coordinates) || ...
            ~ismember(size(Coordinates,2),[2 3]) || ...
            any(~isfinite(Coordinates(:)))
        error('GSPTB:InvalidCoordinates', ...
            ['Coordinates must be an N-by-2 or N-by-3 finite ' ...
             'numeric array of real numbers!']);
    end

    % Checks symmetry of the Adjacency
    Tolerance = 1e-12 * max(1,max(abs(Adjacency(:))));

    if max(abs(Adjacency-Adjacency.'),[],'all') > Tolerance
        error('GSPTB:InvalidAdjacency', ...
            'Adjacency must be symmetric!');
    end

    % Enforces non-negative weights
    if any(Adjacency(:) < 0)
        error('GSPTB:InvalidAdjacency', ...
            'Adjacency must contain non-negative edge weights!');
    end

    % Gets only the upper right triangle
    UpperAdjacency = triu(Adjacency,1);
    [Source,Target,Weights] = find(UpperAdjacency);

    % Validations for NodeColors
    if ~isnumeric(NodeColors) || isempty(NodeColors) || ~isreal(NodeColors)
        error('GSPTB:InvalidNodeColor', ...
            'NodeColors must be a non-empty numeric array of real numbers!');
    end
    
    if isequal(size(NodeColors),[1 3])
    
        if any(~isfinite(NodeColors)) || ...
                any(NodeColors < 0 | NodeColors > 1)
            error('GSPTB:InvalidNodeColor', ...
                'RGB values in NodeColors must lie between 0 and 1!');
        end
    
    elseif iscolumn(NodeColors) && numel(NodeColors) == N
    
        if any(~isfinite(NodeColors))
            error('GSPTB:InvalidNodeColor', ...
                'Scalar node-color values must be finite!');
        end
    
    elseif isequal(size(NodeColors),[N 3])
    
        if any(~isfinite(NodeColors(:))) || ...
                any(NodeColors(:) < 0 | NodeColors(:) > 1)
            error('GSPTB:InvalidNodeColor', ...
                ['Per-node RGB values must be finite and lie ' ...
                 'between 0 and 1!']);
        end
    
    else
        error('GSPTB:InvalidNodeColor', ...
            ['NodeColors must be a 1-by-3 RGB triplet, an N-by-1 ' ...
             'vector of scalar color values, or an N-by-3 RGB array!']);
    end

    % Validations for NodeSizes
    if ~isnumeric(NodeSizes) || isempty(NodeSizes) || ...
        any(~isfinite(NodeSizes(:))) || ~isreal(NodeSizes) || ...
        any(NodeSizes(:) <= 0)
        error('GSPTB:InvalidNodeSize', ...
            'NodeSizes must contain positive finite real numeric values!');
    end
    
    if ~isscalar(NodeSizes) && ...
            ~(iscolumn(NodeSizes) && numel(NodeSizes) == N)
        error('GSPTB:DimensionMismatch', ...
            ['NodeSizes must be scalar or an N-by-1 vector ' ...
             'containing one size per node!']);
    end
    
    NodeSizes = NodeSizes(:);

    % Validations for EdgeColor
    if ~isnumeric(EdgeColor) || ...
        ~isequal(size(EdgeColor),[1 3]) || ~isreal(EdgeColor) || ...
        any(~isfinite(EdgeColor)) || ...
        any(EdgeColor < 0 | EdgeColor > 1)
        error('GSPTB:InvalidEdgeColor', ...
            ['EdgeColor must be a finite 1-by-3 RGB triplet ' ...
             'with real values between 0 and 1!']);
    end

    % Validations for EdgeWidthRange
    if ~isnumeric(EdgeWidthRange) || ...
        ~isequal(size(EdgeWidthRange),[1 2]) || ...
        any(~isfinite(EdgeWidthRange)) || ~isreal(EdgeWidthRange) || ...
        any(EdgeWidthRange <= 0) || ...
        EdgeWidthRange(1) > EdgeWidthRange(2)
        error('GSPTB:InvalidEdgeWidth', ...
            ['EdgeWidthRange must contain two positive finite ' ...
             'real values in non-decreasing order!']);
    end

    % Extracts the widths for plotting graph edges
    if isempty(Weights)
        Widths = zeros(0,1);
    elseif max(Weights) == min(Weights)
        Widths = repmat(mean(EdgeWidthRange),size(Weights));
    else
        Widths = EdgeWidthRange(1) + ...
            (Weights-min(Weights)) ./ ...
            (max(Weights)-min(Weights)) .* ...
            diff(EdgeWidthRange);
    end

    % Keeps track of axis hold status pre-plotting
    WasHeld = ishold(ax);
    hold(ax,'on');
    
    % Number of edges to plot
    E = numel(Weights);

    % Initializes the object that will contain all edges
    hEdges = gobjects(E,1);
    
    % Same process for all edges
    for e = 1:E

        % Gets nodes' location
        Nodes = [Source(e),Target(e)];
    
        if size(Coordinates,2) == 2
            hEdges(e) = line(ax, ...
                Coordinates(Nodes,1), ...
                Coordinates(Nodes,2), ...
                'Color',EdgeColor, ...
                'LineWidth',Widths(e));
        else
            hEdges(e) = line(ax, ...
                Coordinates(Nodes,1), ...
                Coordinates(Nodes,2), ...
                Coordinates(Nodes,3), ...
                'Color',EdgeColor, ...
                'LineWidth',Widths(e));
        end
    end
    
    if size(Coordinates,2) == 2
        hNodes = scatter(ax, ...
            Coordinates(:,1),Coordinates(:,2), ...
            NodeSizes,NodeColors,'filled');
    else
        hNodes = scatter3(ax, ...
            Coordinates(:,1),Coordinates(:,2),Coordinates(:,3), ...
            NodeSizes,NodeColors,'filled');
    end
    
    if ~WasHeld
        hold(ax,'off');
    end
    
    axis(ax,'equal');
    axis(ax,'off');

    h.Figure = f;
    h.Axes = ax;
    h.Edges = hEdges;
    h.Nodes = hNodes;
    h.Source = Source;
    h.Target = Target;
    h.EdgeWeights = Weights;
end