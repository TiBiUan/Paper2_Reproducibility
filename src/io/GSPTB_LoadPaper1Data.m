function [Data] = GSPTB_LoadPaper1Data(Cfg,Root)

    Data.Parents = loadVariable(fullfile(Root,'Parents.mat'),'Parents');
    Data.Parents = Data.Parents(Cfg.Data.ToLoad,:);

    Data.SubjectIndices = loadVariable(fullfile(Root,'Subject_Indices.mat'),'Subject_Indices');
    Data.SubjectIndices = Data.SubjectIndices(Cfg.Data.ToLoad);

    Data.RegionDistances = loadVariable(fullfile(Root,'d_scale3.mat'),'d_scale3');

    tmp = loadVariable(fullfile(Root,'Factor_Scores.mat'),'Factor_Scores');
    Data.FactorScores.MentalHealth = tmp(Cfg.Data.ToLoad,1);
    Data.FactorScores.Cognition = tmp(Cfg.Data.ToLoad,2);
    Data.FactorScores.ProcessingSpeed = tmp(Cfg.Data.ToLoad,3);
    Data.FactorScores.SubstanceUse = tmp(Cfg.Data.ToLoad,4);
    
    Data.Adj = loadVariable(fullfile(Root,'Adj.mat'),'Adj');
    Data.Adj = Data.Adj(Cfg.Data.ToLoad);

    Data.DegMat = loadVariable(fullfile(Root,'DegMat.mat'),'DegMat');
    Data.DegMat = Data.DegMat(Cfg.Data.ToLoad);
    
    Data.U = loadVariable(fullfile(Root,'U.mat'),'U');
    Data.U = Data.U(Cfg.Data.ToLoad);
   
    Data.Lambda = loadVariable(fullfile(Root,'Lambda.mat'),'Lambda');
    Data.Lambda = Data.Lambda(Cfg.Data.ToLoad);

    Data.NumSubjects = numel(Data.SubjectIndices);
    Data.NumRegions = Cfg.Atlas.NumRegions;

    Data.FamilyKernel = GSPTB_IsSameFamily(Data.Parents);

    Data.FD_TimeCourses = loadVariable(fullfile(Root,'FD_tasks.mat'),'FD_tasks');
    Data.FD_TimeCourses = cellfun(@(x) x(Cfg.Data.ToLoad),Data.FD_TimeCourses,'UniformOutput',false);

    Data.MeanFD = loadVariable(fullfile(Root,'MFD.mat'),'MFD');
    Data.MeanFD = Data.MeanFD(Cfg.Data.ToLoad,:);

    Data.MeanFD_NonScrubbedOnly = loadVariable(fullfile(Root,'MFD_NS.mat'),'MFD_NS');
    Data.MeanFD_NonScrubbedOnly = Data.MeanFD_NonScrubbedOnly(Cfg.Data.ToLoad,:);

    Data.PercentageScrubbed = loadVariable(fullfile(Root,'P_SC.mat'),'P_SC');
    Data.PercentageScrubbed = Data.PercentageScrubbed(Cfg.Data.ToLoad,:);

    validatePaper1Data(Data,Cfg);
end

function Value = loadVariable(FilePath,VariableName)

    if ~isfile(FilePath)
        error('GSPTB:MissingFile', ...
            'Required file not found: %s',FilePath);
    end

    Loaded = load(FilePath,VariableName);

    if ~isfield(Loaded,VariableName)
        error('GSPTB:MissingVariable', ...
            'Variable "%s" was not found in %s.', ...
            VariableName,FilePath);
    end

    Value = Loaded.(VariableName);
end

function validatePaper1Data(Data,Cfg)

    S = Data.NumSubjects;
    R = Cfg.Atlas.NumRegions;

    if size(Data.Parents,1) ~= S
        error('GSPTB:SubjectMismatch', ...
            'Parents and SubjectIndices contain different subject counts.');
    end

    if size(Data.FactorScores.MentalHealth,1) ~= S || size(Data.FactorScores.Cognition,1) ~= S ||...
            size(Data.FactorScores.ProcessingSpeed,1) ~= S || size(Data.FactorScores.SubstanceUse,1) ~= S

        error('GSPTB:SubjectMismatch', ...
            'Factor scores must contain one entry per subject.');
    end

    if numel(Data.Adj) ~= S || numel(Data.U) ~= S || ...
            numel(Data.Lambda) ~= S
        error('GSPTB:SubjectMismatch', ...
            'Structural data do not match the number of subjects.');
    end

    if ~isequal(size(Data.RegionDistances),[R R])
        error('GSPTB:RegionMismatch', ...
            'The region-distance matrix must be %d-by-%d.',R,R);
    end

    for s = 1:S
        if ~isequal(size(Data.Adj{s}),[R R])
            error('GSPTB:RegionMismatch', ...
                'Adj{%d} is not %d-by-%d.',s,R,R);
        end

        if ~isequal(size(Data.U{s}),[R R])
            error('GSPTB:RegionMismatch', ...
                'U{%d} is not %d-by-%d.',s,R,R);
        end

        if numel(Data.Lambda{s}) ~= R
            error('GSPTB:RegionMismatch', ...
                'Lambda{%d} does not contain %d eigenvalues.',s,R);
        end
    end
end