function [Results] = GSPTB_RankSum_TaskContrasts(Metric,Cfg)

    R = size(Metric{1},1);
    S = size(Metric{1},2);

    % We want to test each harmonic in terms of potential differences across
    % tasks
    for c = 1:R

        % Creates the data table
        ws_oi = [];
        
        for t = 1:length(Cfg.Tasks.Names)
            ws_oi = [ws_oi;Metric{t}(c,:)'];
        end
        
        task_oi = repmat(1:length(Cfg.Tasks.Names),S,1);
        task_oiv = task_oi(:);
        
        % Post-hoc testing between pairs of tasks
        for t1 = 1:length(Cfg.Tasks.Names)
            for t2 = 1:length(Cfg.Tasks.Names)
                if t1 ~= t2
                    [Results.pValue(c,t1,t2),~,stats] = ranksum(ws_oi(task_oiv==t1),ws_oi(task_oiv==t2));
                    Results.sign(c,t1,t2) = sign(stats.zval);
                    Results.zStat(c,t1,t2) = stats.zval;
                end
            end
        end
    end
end