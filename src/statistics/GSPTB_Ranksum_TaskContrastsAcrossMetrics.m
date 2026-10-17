function [TaskContrastResults_RankSum] = GSPTB_Ranksum_TaskContrastsAcrossMetrics(DynMet,Cfg)

    MeanIntensity = cellfun(@(x) x.MeanIntensity,DynMet,'UniformOutput',false);
    TaskContrastResults_RankSum.MeanIntensity = GSPTB_RankSum_TaskContrasts(MeanIntensity,Cfg);

    EntryRate = cellfun(@(x) x.EntryRate,DynMet,'UniformOutput',false);
    TaskContrastResults_RankSum.EntryRate = GSPTB_RankSum_TaskContrasts(EntryRate,Cfg);

    Duration = cellfun(@(x) x.Duration,DynMet,'UniformOutput',false);
    TaskContrastResults_RankSum.Duration = GSPTB_RankSum_TaskContrasts(Duration,Cfg);
end