function [DynamicMetrics,SparsePredictionResults] = GSPTB_RunCohortDynamicAnalysis(Cfg,Paths,U_R,task_index)

    % Loads the data
    X_OI = GSPTB_LoadTaskData(Paths.Functional,Cfg.Tasks.Names{task_index},...
        Cfg.Atlas.CorticalIndices,Cfg.Atlas.IsCortexOnly);

    % Number of subjects
    S = length(U_R);
        
    SparsePredictionResults = cell(1,S);

    % Analyses per se: (1) sparse representation, (2) dynamic metrics
    for s = 1:S

        disp(['Dynamic analysis for subject ',num2str(s),'...']);
        
        % We perform the prediction for the current task/subject data
        SparsePredictionResults{s} = GSPTB_RunSparseActivityPrediction(X_OI{s}{1}',X_OI{s}{2}',U_R{s});
        
        % Global weights used to generate the metrics
        DynamicMetrics.MeanIntensity(:,s) = SparsePredictionResults{s}.MeanIntensity;
        [DynamicMetrics.EntryRate(:,s),DynamicMetrics.Duration(:,s)] = ...
            GSPTB_ExtractDynamicMetrics(SparsePredictionResults{s}.BinWeights,Cfg.Data.TR);
        
        % Positive-only and negative-only weights used instead
        DynamicMetrics.PosOnly.MeanIntensity(:,s) = SparsePredictionResults{s}.PosOnly.MeanIntensity;
        DynamicMetrics.NegOnly.MeanIntensity(:,s) = SparsePredictionResults{s}.NegOnly.MeanIntensity;
        
        [DynamicMetrics.PosOnly.EntryRate(:,s),DynamicMetrics.PosOnly.Duration(:,s)] = ...
            GSPTB_ExtractDynamicMetrics(SparsePredictionResults{s}.PosOnly.BinWeights,Cfg.Data.TR);
        
        [DynamicMetrics.NegOnly.EntryRate(:,s),DynamicMetrics.NegOnly.Duration(:,s)] = ...
            GSPTB_ExtractDynamicMetrics(SparsePredictionResults{s}.NegOnly.BinWeights,Cfg.Data.TR);
    end









end