function [h] = GSPTB_ExportFigure(FigureHandle,FileName,Format,Resolution,BackgroundColor)

    % Validate figure
    if isempty(FigureHandle) || ...
            ~isgraphics(FigureHandle,'figure')
        error('GSPTB:InvalidFigure', ...
            'FigureHandle must be a valid figure handle!');
    end

    % Validate filename
    if ischar(FileName)
        FileName = string(FileName);
    end

    if ~isstring(FileName) || ~isscalar(FileName) || ...
            ismissing(FileName) || strlength(FileName) == 0
        error('GSPTB:InvalidFileName', ...
            'FileName must be a non-empty text scalar!');
    end

    % Validate format
    Format = lower(string(Format));
    ValidFormats = ["png","pdf","svg","eps","jpg","jpeg","tif","tiff"];

    if ~isscalar(Format) || ~ismember(Format,ValidFormats)
        error('GSPTB:InvalidFormat', ...
            ['Format must be ''png'', ''pdf'', ''svg'', ''eps'', ' ...
             '''jpg'', or ''tif''!']);
    end

    % Validate resolution
    if ~isnumeric(Resolution) || ~isreal(Resolution) || ...
            ~isscalar(Resolution) || ...
            ~isfinite(Resolution) || Resolution <= 0
        error('GSPTB:InvalidResolution', ...
            'Resolution must be a positive finite scalar!');
    end

    % Validate background color
    if ischar(BackgroundColor)
        BackgroundColor = string(BackgroundColor);
    end

    if isstring(BackgroundColor)
        if ~isscalar(BackgroundColor) || ...
                lower(BackgroundColor) ~= "none"
            error('GSPTB:InvalidBackgroundColor', ...
                ['BackgroundColor must be ''none'' or a finite ' ...
                 '1-by-3 RGB triplet!']);
        end
        BackgroundColor = "none";

    elseif ~isnumeric(BackgroundColor) || ...
            ~isreal(BackgroundColor) || ...
            ~isequal(size(BackgroundColor),[1 3]) || ...
            any(~isfinite(BackgroundColor)) || ...
            any(BackgroundColor < 0 | BackgroundColor > 1)

        error('GSPTB:InvalidBackgroundColor', ...
            ['BackgroundColor must be ''none'' or a finite ' ...
             '1-by-3 RGB triplet!']);
    end

    % Normalize common format aliases
    if Format == "jpeg"
        Format = "jpg";
    elseif Format == "tiff"
        Format = "tif";
    end

    % Remove an existing extension so the requested format controls output
    [Folder,BaseName,~] = fileparts(FileName);

    if strlength(Folder) == 0
        OutputFile = BaseName + "." + Format;
    else
        OutputFile = fullfile(Folder,BaseName + "." + Format);
    end

    % Create destination folder when needed
    OutputFolder = fileparts(OutputFile);

    if strlength(OutputFolder) > 0 && ~isfolder(OutputFolder)
        mkdir(OutputFolder);
    end

    % Raster and vector formats require slightly different handling
    RasterFormats = ["png","jpg","tif"];

    if ismember(Format,RasterFormats)
        exportgraphics(FigureHandle,OutputFile, ...
            'Resolution',Resolution, ...
            'BackgroundColor',BackgroundColor);
    else
        exportgraphics(FigureHandle,OutputFile, ...
            'ContentType','vector', ...
            'BackgroundColor',BackgroundColor);
    end

    % Return export metadata
    h.Figure = FigureHandle;
    h.FileName = OutputFile;
    h.Format = Format;
    h.Resolution = Resolution;
    h.BackgroundColor = BackgroundColor;
end