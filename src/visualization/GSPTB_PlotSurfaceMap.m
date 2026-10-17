function [h] = GSPTB_PlotSurfaceMap(ax,Values,AtlasFile,SurfaceFile,...
    ColorMap,ColorLimits,View)

    % GSPTB_PlotSurfaceMap
    %
    % Maps regional values onto a labelled cortical atlas and renders the
    % resulting map on a cortical surface using BrainNet Viewer.
    %
    % INPUTS
    %   ax
    %       Destination axes. Pass [] to create a new figure and axes.
    %
    %   Values
    %       N-by-1 vector containing one value per cortical region.
    %
    %   AtlasFile
    %       Cortex-only labelled NIfTI atlas. Its nonzero labels must be
    %       the consecutive integers 1:N.
    %
    %   SurfaceFile
    %       BrainNet Viewer cortical-surface file (.nv).
    %
    %   ColorMap
    %       K-by-3 RGB colormap with values between 0 and 1.
    %
    %   ColorLimits
    %       1-by-2 vector containing the lower and upper color limits.
    %
    %   View
    %       One of: 'left', 'right', 'front', 'back', 'top', 'bottom',
    %       or a 1-by-2 numeric vector containing [azimuth elevation].
    %
    % OUTPUT
    %   h
    %       Structure containing the destination figure, axes, copied
    %       graphics objects, regional values, and input-file metadata.

    Views.left.Az    = -90;
    Views.left.El    =   0;
    Views.left.Up    = [0 0 1];
    
    Views.right.Az   =  90;
    Views.right.El   =   0;
    Views.right.Up   = [0 0 1];
    
    Views.front.Az   = 180;
    Views.front.El   =   0;
    Views.front.Up   = [0 0 1];
    
    Views.back.Az    =   360;
    Views.back.El    =   0;
    Views.back.Up    = [0 0 1];
    
    Views.top.Az     = 270;
    Views.top.El     =  90;
    Views.top.Up     = [1 0 0];
    
    Views.bottom.Az  = -90;
    Views.bottom.El  = -90;
    Views.bottom.Up  = [1 0 0];


    %% Create or validate destination axes

    if isempty(ax)

        f = figure;
        ax = axes(f);

    elseif ~isscalar(ax) || ~isgraphics(ax,'axes')
        error('GSPTB:InvalidAxes', ...
            'ax must be empty or a valid scalar axes handle!');
    else
        f = ancestor(ax,'figure');
    end


    %% Validate regional values

    if isempty(Values) || ~isnumeric(Values) || ...
            ~isreal(Values) || ~isvector(Values) || ...
            any(~isfinite(Values(:)))

        error('GSPTB:InvalidSurfaceValues', ...
            ['Values must be a non-empty finite real ' ...
             'numeric vector!']);
    end

    Values = Values(:);
    N = numel(Values);


    %% Validate files

    AtlasFile = localValidateFile( ...
        AtlasFile, ...
        'GSPTB:MissingAtlas', ...
        'AtlasFile must identify an existing NIfTI atlas file!');

    SurfaceFile = localValidateFile( ...
        SurfaceFile, ...
        'GSPTB:MissingSurface', ...
        'SurfaceFile must identify an existing BrainNet .nv file!');


    %% Validate colormap

    if ~isnumeric(ColorMap) || ~isreal(ColorMap) || ...
            ~ismatrix(ColorMap) || size(ColorMap,2) ~= 3 || ...
            size(ColorMap,1) < 2 || ...
            any(~isfinite(ColorMap(:))) || ...
            any(ColorMap(:) < 0 | ColorMap(:) > 1)

        error('GSPTB:InvalidColorMap', ...
            ['ColorMap must be a finite K-by-3 RGB array ' ...
             'with values between 0 and 1!']);
    end


    %% Validate color limits

    if ~isnumeric(ColorLimits) || ~isreal(ColorLimits) || ...
            ~isequal(size(ColorLimits),[1 2]) || ...
            any(~isfinite(ColorLimits)) || ...
            ColorLimits(1) >= ColorLimits(2)

        error('GSPTB:InvalidColorLimits', ...
            ['ColorLimits must contain two finite values in ' ...
             'strictly increasing order!']);
    end


    %% Validate view

    if ischar(View)
        View = string(View);
    end

    IsNamedView = isstring(View) && isscalar(View) && ...
        ~ismissing(View) && strlength(View) > 0;

    IsNumericView = isnumeric(View) && isreal(View) && ...
        isequal(size(View),[1 2]) && all(isfinite(View));

    if ~IsNamedView && ~IsNumericView

        error('GSPTB:InvalidView', ...
            ['View must be a named view or a finite 1-by-2 ' ...
             'azimuth/elevation vector!']);
    end

    if IsNamedView

        View = lower(View);

        ValidViews = [
            "left"
            "right"
            "front"
            "back"
            "top"
            "bottom"
        ];

        if ~ismember(View,ValidViews)

            error('GSPTB:InvalidView', ...
                ['Named View must be left, right, front, back, ' ...
                 'top, or bottom!']);

        end
    end


    %% Verify dependencies

    if exist('BrainNet_MapCfg','file') ~= 2

        error('GSPTB:MissingDependency', ...
            ['BrainNet Viewer could not be found on the MATLAB path. ' ...
             'Run GSPTB_SetupDependencies before plotting!']);
    end

    if exist('spm_vol','file') ~= 2 || ...
            exist('spm_read_vols','file') ~= 2 || ...
            exist('spm_write_vol','file') ~= 2 || ...
            exist('spm_type','file') ~= 2

        error('GSPTB:MissingDependency', ...
            ['SPM could not be found on the MATLAB path. ' ...
             'Surface plotting requires SPM volume I/O functions!']);
    end


    %% Load atlas

    try

        AtlasHeader = spm_vol(char(AtlasFile));

        if numel(AtlasHeader) ~= 1

            error('GSPTB:InvalidAtlas', ...
                'AtlasFile must contain one three-dimensional volume!');
        end

        Atlas = spm_read_vols(AtlasHeader);

    catch ME

        if strcmp(ME.identifier,'GSPTB:InvalidAtlas')
            rethrow(ME);
        end

        NewException = MException( ...
            'GSPTB:InvalidAtlas', ...
            'The atlas NIfTI file could not be read with SPM.');

        NewException = addCause(NewException,ME);
        throw(NewException);

    end


    %% Validate atlas labels

    if ~isnumeric(Atlas) || ~isreal(Atlas) || ...
            any(~isfinite(Atlas(:))) || any(Atlas(:) < 0)

        error('GSPTB:InvalidAtlas', ...
            ['The atlas must contain finite, non-negative, ' ...
             'real numeric labels!']);
    end

    RoundedAtlas = round(Atlas);

    if any(abs(Atlas(:)-RoundedAtlas(:)) > 1e-8)

        error('GSPTB:InvalidAtlas', ...
            'The atlas must contain integer-valued region labels!');
    end

    Atlas = RoundedAtlas;

    AtlasLabels = unique(Atlas(:));
    AtlasLabels(AtlasLabels == 0) = [];

    if ~isequal(AtlasLabels(:),(1:N).')

        error('GSPTB:AtlasValueMismatch', ...
            ['The atlas must contain exactly the consecutive labels ' ...
             '1:%d to match the %d regional values!'],N,N);
    end


    %% Map regional values into voxel space

    SurfaceVolume = zeros(size(Atlas));

    for r = 1:N
        SurfaceVolume(Atlas == r) = Values(r);
    end


    %% Write temporary NIfTI

    TemporaryMapFile = [tempname '.nii'];

    OutputHeader = AtlasHeader;
    OutputHeader.fname = TemporaryMapFile;
    OutputHeader.dt = [spm_type('float32') 0];
    OutputHeader.pinfo = [1; 0; 0];
    OutputHeader.descrip = 'GSPTB regional surface map';

    try

        spm_write_vol(OutputHeader,SurfaceVolume);

    catch ME

        if isfile(TemporaryMapFile)
            delete(TemporaryMapFile);
        end

        NewException = MException( ...
            'GSPTB:SurfaceMapWriteFailure', ...
            ['The temporary surface-map NIfTI could not be ' ...
             'written with SPM.']);

        NewException = addCause(NewException,ME);
        throw(NewException);
    end


    %% Render in BrainNet's temporary figure

    FiguresBefore = findall(groot,'Type','figure');

    try
        BrainNet_MapCfg( ...
            char(SurfaceFile), ...
            char(TemporaryMapFile));

    catch ME

        if isfile(TemporaryMapFile)
            delete(TemporaryMapFile);
        end

        NewException = MException( ...
            'GSPTB:BrainNetRenderingFailure', ...
            'BrainNet Viewer could not render the surface map.');

        NewException = addCause(NewException,ME);
        throw(NewException);
    end

    FiguresAfter = findall(groot,'Type','figure');
    NewFigures = setdiff(FiguresAfter,FiguresBefore);

    if isempty(NewFigures)
        TemporaryFigure = gcf;
    else
        TemporaryFigure = NewFigures(1);
    end


    %% Identify the BrainNet surface axes

    TemporaryAxes = findall(TemporaryFigure,'Type','axes');

    PatchCounts = arrayfun(@(CurrentAxes) ...
        numel(findall(CurrentAxes,'Type','patch')), ...
        TemporaryAxes);

    [LargestPatchCount,SurfaceAxesIndex] = max(PatchCounts);

    if LargestPatchCount == 0

        if isgraphics(TemporaryFigure)
            close(TemporaryFigure);
        end

        if isfile(TemporaryMapFile)
            delete(TemporaryMapFile);
        end

        error('GSPTB:SurfaceAxesNotFound', ...
            'The rendered cortical-surface axes could not be identified!');
    end

    SourceAxes = TemporaryAxes(SurfaceAxesIndex);


    %% Copy the rendered surface into the destination axes

    WasHeld = ishold(ax);
    hold(ax,'on');

    CopiedObjects = copyobj(allchild(SourceAxes),ax);

    if ~WasHeld
        hold(ax,'off');
    end


    %% Copy relevant axes properties

    ax.XLim = SourceAxes.XLim;
    ax.YLim = SourceAxes.YLim;
    ax.ZLim = SourceAxes.ZLim;

    ax.DataAspectRatio = SourceAxes.DataAspectRatio;
    ax.PlotBoxAspectRatio = SourceAxes.PlotBoxAspectRatio;

    ax.CameraPosition = SourceAxes.CameraPosition;
    ax.CameraTarget = SourceAxes.CameraTarget;
    ax.CameraUpVector = SourceAxes.CameraUpVector;
    ax.CameraViewAngle = SourceAxes.CameraViewAngle;

    ax.CLim = ColorLimits;

    colormap(ax,ColorMap);

    axis(ax,'off');
    axis(ax,'vis3d');


    %% Makes sure right is on the right for "front" view

    if IsNamedView && View == "front"

        SurfacePatches = findall(ax,'Type','patch');
    
        for p = 1:numel(SurfacePatches)
    
            Vertices = SurfacePatches(p).Vertices;
    
            if isempty(Vertices)
                continue;
            end
    
            % Reflect around the centre of the surface's left-right extent.
            ReflectionCentre = ...
                (min(Vertices(:,1)) + max(Vertices(:,1))) / 2;
    
            Vertices(:,1) = ...
                2*ReflectionCentre - Vertices(:,1);
    
            SurfacePatches(p).Vertices = Vertices;
    
            % Makes lighting robust to the reflected surface orientation.
            SurfacePatches(p).BackFaceLighting = 'reverselit';
    
        end
    end


    %% Apply requested view

    ZoomFactor = 1.15;
    V = Views.(lower(View));

    view(ax,V.Az,V.El);
    camup(ax,V.Up);

    if IsNamedView && View == "bottom"
        ax.XDir = 'reverse';
    else
        ax.XDir = 'normal';
    end

    axis(ax,'vis3d');
    axis(ax,'tight');

    camzoom(ax,ZoomFactor);

    % Ensures that camera properties are updated before lighting.
    drawnow;


    %% Recreate view-consistent lighting

    delete(findall(ax,'Type','light'));

    % Direction from the camera target toward the camera. This creates
    % headlight-like illumination without relying on camlight, which can
    % fail when MATLAB temporarily reports a non-finite camera position.
    LightDirection = ax.CameraPosition - ax.CameraTarget;

    if any(~isfinite(LightDirection)) || norm(LightDirection) == 0

        % Defensive fallback. This should rarely be required.
        LightDirection = [0 0 1];
    else
        LightDirection = LightDirection ./ norm(LightDirection);
    end

    hLight = light(ax, ...
        'Style','infinite', ...
        'Position',LightDirection);

    lighting(ax,'gouraud');


    %% Remove temporary files and BrainNet figure

    if isgraphics(TemporaryFigure)
        close(TemporaryFigure);
    end

    if isfile(TemporaryMapFile)
        delete(TemporaryMapFile);
    end


    %% Return handles and metadata

    h.Figure = f;
    h.Axes = ax;
    h.Objects = CopiedObjects;
    h.Patches = findall(ax,'Type','patch');
    h.Light = hLight;
    h.Values = Values;
    h.AtlasFile = AtlasFile;
    h.SurfaceFile = SurfaceFile;
    h.ColorMap = ColorMap;
    h.ColorLimits = ColorLimits;
    h.View = View;

end


function FileName = localValidateFile( ...
    FileName,ErrorIdentifier,ErrorMessage)

    if ischar(FileName)
        FileName = string(FileName);
    end

    if ~isstring(FileName) || ~isscalar(FileName) || ...
            ismissing(FileName) || strlength(FileName) == 0

        error(ErrorIdentifier,'%s',ErrorMessage);
    end

    if ~isfile(FileName)

        error(ErrorIdentifier, ...
            '%s Received: %s',ErrorMessage,FileName);
    end
end