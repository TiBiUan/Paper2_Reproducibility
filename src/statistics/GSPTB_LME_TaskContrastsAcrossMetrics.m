function [TaskContrasts] = GSPTB_LME_TaskContrastsAcrossMetrics(Stats,Cfg)

    TaskContrasts.MeanIntensity = GSPTB_LME_TaskContrasts(Stats.MeanIntensity.MixedEffects.Model2.Model,...
        Cfg.TaskContrasts.TaskFactorName,Cfg.Statistics.TaskContrasts.AlphaLevel);

    TaskContrasts.EntryRate = GSPTB_LME_TaskContrasts(Stats.EntryRate.MixedEffects.Model2.Model,...
        Cfg.TaskContrasts.TaskFactorName,Cfg.Statistics.TaskContrasts.AlphaLevel);
    
    TaskContrasts.Duration = GSPTB_LME_TaskContrasts(Stats.Duration.MixedEffects.Model2.Model,...
        Cfg.TaskContrasts.TaskFactorName,Cfg.Statistics.TaskContrasts.AlphaLevel);
end