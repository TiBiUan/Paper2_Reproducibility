% X has size R x 1, R number of brain regions/harmonics and should be
% spectral index (or another monotonically increasing quantity)
% Y has size N x R, with N the number of data points
% ErrorType specifies with type of error to display between "STD", "SEM",
% "CI" or "NONE"
% Color specifies the color of the curve (1 x 3 vector)
% Alpha specifies the transparency of the error surface
function [h] = GSPTB_PlotSpectralCurve(ax,X,Y,ErrorType,Color,Alpha)
    
    % If we provide an empty axis handle, then we create what we need
    if isempty(ax)
        f = figure;
        ax = axes(f);

    % Else, we only get the underlying figure handle from what was given as
    % input axis
    else
        f = ancestor(ax,'figure');
    end

    % Checks for empty inputs
    if isempty(X) || isempty(Y)
        error('GSPTB:EmptyInput','X and Y must be non-empty!');
    end

    % Checks dimensionality of the inputs
    if ~isvector(X)
        error('GSPTB:InvalidX','X must be a vector!');
    end
    
    X = X(:).';
    
    if size(Y,2) ~= numel(X)
        error('GSPTB:DimensionMismatch', ...
            'Y has %d columns, but X contains %d elements!',size(Y,2),numel(X));
    end

    % Mean curve
    MU  = mean(Y,1,'omitnan');

    % Uncertainties, accounting for possible NaN values
    STD = std(Y,0,1,'omitnan');
    N = sum(~isnan(Y),1);
    SEM = STD./sqrt(N);
    CI = tinv(0.975,N-1) .* SEM;
    CI(N < 2) = NaN;

    % Preserve the initial hold state
    WasHeld = ishold(ax);
    hold(ax,'on');
    
    % Builds the uncertainty patch
    switch upper(ErrorType)
        case 'STD'
            Error = STD;
        case 'SEM'
            Error = SEM;
        case 'CI'
            Error = CI;
        case 'NONE'
            Error = [];
        otherwise
            error('GSPTB:InvalidErrorType', ...
                'ErrorType must be ''STD'', ''SEM'', ''CI'', or ''None''!');
    end

    % Plots the uncertainty surfaces if required
    if isempty(Error)
        hPatch = gobjects(0);
    else
        hPatch = fill(ax,[X fliplr(X)],[MU+Error fliplr(MU-Error)],Color,...
            'EdgeColor','None','FaceAlpha',Alpha);
    end

    % Plot the mean
    hLine = plot(ax,X,MU,'Color',Color,'LineWidth',2,'LineStyle','-');

    % Restore the original hold state
    if ~WasHeld
        hold(ax,'off');
    end

    % Aesthetics
    xlim(ax,[min(X),max(X)]);
    set(ax,'Box','off');

    % Prepares object to return
    h.Figure = f;
    h.Axes = ax;
    h.Uncertainty = hPatch;
    h.Line = hLine;
end