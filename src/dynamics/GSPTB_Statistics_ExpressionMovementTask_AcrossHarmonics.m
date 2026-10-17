function [CorrMat,Results] = GSPTB_Statistics_ExpressionMovementTask_AcrossHarmonics(Motion,Metric)

    % Number of tasks
    n_tasks = size(Motion,2);

    % Number of subjects
    S = size(Motion,1);

    % Number of harmonics
    R = size(Metric{1},1);

    % We first compute the matrix of movement/metric correlations across
    % harmonics and tasks
    CorrMat = NaN(R,n_tasks);

    for t = 1:n_tasks    
        for c = 1:R
            CorrMat(c,t) = corr(Motion(:,t),Metric{t}(c,:)');
        end
    end

    % Next, we perform mixed modelling with three models: only including
    % task, including task and movement, and also their interaction, to
    % describe the metric on top of a random subject intercept
    tmp_task = repmat(1:n_tasks,S,1);
    tmp_subj = repmat((1:S)',1,n_tasks);

    Subject_vec = categorical(tmp_subj(:));
    Task_vec = categorical(tmp_task(:));
    Motion_vec = Motion(:);
    Motion_centered = Motion_vec - mean(Motion_vec);

    DataTab = table(Subject_vec,Task_vec,Motion_centered,...
        zeros(numel(Motion_vec),1),'VariableNames',...
    {'Subject','Task','Motion','Metric'});

    % Same process for each harmonic...
    for c = 1:R

        c

        % Adds the relevant Metric data in the table
        MetData = NaN(S,n_tasks);

        for t = 1:n_tasks
            MetData(:,t) = Metric{t}(c,:)';
        end

        Metric_vec = MetData(:);
        DataTab.Metric = Metric_vec;

        % First model
        lme_0 = fitlme(DataTab,'Metric ~ Task + (1|Subject)',...
        'FitMethod','ML');

        % Second model
        lme_1 = fitlme(DataTab,'Metric ~ Task + Motion + (1|Subject)',...
        'FitMethod','ML');
        
        % Third model
        lme_2 = fitlme(DataTab,'Metric ~ Task * Motion + (1|Subject)',...
        'FitMethod','ML');

        % Associated approach to derive an omnibus statistic for the task
        % contrast
        Stats0 = anova(lme_0,'DFMethod','satterthwaite');
        Stats1 = anova(lme_1,'DFMethod','satterthwaite');
        Stats2 = anova(lme_2,'DFMethod','satterthwaite');

        % Comparison between both models
        tmp_comp_01 = compare(lme_0,lme_1);
        tmp_comp_12 = compare(lme_1,lme_2);

        % Now sampling the resulting outputs to return them
        idx_T0 = strcmp(string(Stats0.Term),'Task');

        idx_T1 = strcmp(string(Stats1.Term),'Task');
        idx_M1 = strcmp(lme_1.Coefficients.Name,'Motion');

        idx_T2 = strcmp(string(Stats2.Term),'Task');
        idx_M2 = strcmp(lme_2.Coefficients.Name,'Motion');
        idx_TM2 = strcmp(string(Stats2.Term),'Task:Motion');

        % Results for the first model (minimal)
        Results.Model1.Task.OmnibusF(c)   = Stats0.FStat(idx_T0);
        Results.Model1.Task.DF1(c) = Stats0.DF1(idx_T0);
        Results.Model1.Task.DF2(c) = Stats0.DF2(idx_T0);
        Results.Model1.Task.OmnibusP(c)   = Stats0.pValue(idx_T0);

        Results.Model1.LogLikelihood(c) = lme_0.LogLikelihood;
        Results.Model1.AIC(c) = lme_0.ModelCriterion{1,1};
        Results.Model1.BIC(c) = lme_0.ModelCriterion{1,2};

        % Results for the second model (task and movement)
        Results.Model2.Model{c} = lme_1;

        Results.Model2.Task.OmnibusF(c)   = Stats1.FStat(idx_T1);
        Results.Model2.Task.DF1(c) = Stats1.DF1(idx_T1);
        Results.Model2.Task.DF2(c) = Stats1.DF2(idx_T1);
        Results.Model2.Task.OmnibusP(c)   = Stats1.pValue(idx_T1);
        
        Results.Model2.Motion.Slope(c) = lme_1.Coefficients.Estimate(idx_M1);
        Results.Model2.Motion.SE(c) = lme_1.Coefficients.SE(idx_M1);
        Results.Model2.Motion.t(c) = lme_1.Coefficients.tStat(idx_M1);
        Results.Model2.Motion.p(c) = lme_1.Coefficients.pValue(idx_M1);

        [~,~,tmp] = fixedEffects(lme_1,'DFMethod','satterthwaite');
        Results.Model2.Motion.CI95(c,:) = [tmp.Lower(idx_M1), tmp.Upper(idx_M1)];

        Results.Model2.LogLikelihood(c) = lme_1.LogLikelihood;
        Results.Model2.AIC(c) = lme_1.ModelCriterion{1,1};
        Results.Model2.BIC(c) = lme_1.ModelCriterion{1,2};

        % Results for the third model (full)
        Results.Model3.Task.OmnibusF(c)   = Stats2.FStat(idx_T2);
        Results.Model3.Task.DF1(c) = Stats2.DF1(idx_T2);
        Results.Model3.Task.DF2(c) = Stats2.DF2(idx_T2);
        Results.Model3.Task.OmnibusP(c)   = Stats2.pValue(idx_T2);
        
        Results.Model3.Motion.Slope(c) = lme_2.Coefficients.Estimate(idx_M2);
        Results.Model3.Motion.SE(c) = lme_2.Coefficients.SE(idx_M2);
        Results.Model3.Motion.t(c) = lme_2.Coefficients.tStat(idx_M2);
        Results.Model3.Motion.p(c) = lme_2.Coefficients.pValue(idx_M2);

        [~,~,tmp] = fixedEffects(lme_2,'DFMethod','satterthwaite');
        Results.Model3.Motion.CI95(c,:) = [tmp.Lower(idx_M2), tmp.Upper(idx_M2)];

        Results.Model3.TaskxMotion.OmnibusF(c) = Stats2.FStat(idx_TM2);
        Results.Model3.TaskxMotion.DF1(c)      = Stats2.DF1(idx_TM2);
        Results.Model3.TaskxMotion.DF2(c)      = Stats2.DF2(idx_TM2);
        Results.Model3.TaskxMotion.OmnibusP(c) = Stats2.pValue(idx_TM2);

        Results.Model3.LogLikelihood(c) = lme_2.LogLikelihood;
        Results.Model3.AIC(c) = lme_2.ModelCriterion{1,1};
        Results.Model3.BIC(c) = lme_2.ModelCriterion{1,2};
        
        % Comparison between models
        Results.Model1vs2.LRStat(c) = tmp_comp_01{2,6};
        Results.Model1vs2.p(c) = tmp_comp_01{2,8};

        Results.Model2vs3.LRStat(c) = tmp_comp_12{2,6};
        Results.Model2vs3.p(c) = tmp_comp_12{2,8};
    end
end