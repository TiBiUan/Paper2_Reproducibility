function [MovementExpressionStatistics] = GSPTB_RunMovementExpressionAnalysis(Mot,DynMet)
    
    MeanIntensity = cellfun(@(x) x.MeanIntensity,DynMet,'UniformOutput',false);

    [MovementExpressionStatistics.MeanIntensity.Correlation,...
        MovementExpressionStatistics.MeanIntensity.MixedEffects] = ...
        GSPTB_Statistics_ExpressionMovementTask_AcrossHarmonics(Mot,MeanIntensity);

    EntryRate = cellfun(@(x) x.EntryRate,DynMet,'UniformOutput',false);

    [MovementExpressionStatistics.EntryRate.Correlation,...
        MovementExpressionStatistics.EntryRate.MixedEffects] = ...
        GSPTB_Statistics_ExpressionMovementTask_AcrossHarmonics(Mot,EntryRate);

    Duration = cellfun(@(x) x.Duration,DynMet,'UniformOutput',false);

    [MovementExpressionStatistics.Duration.Correlation,...
        MovementExpressionStatistics.Duration.MixedEffects] = ...
        GSPTB_Statistics_ExpressionMovementTask_AcrossHarmonics(Mot,Duration);
end