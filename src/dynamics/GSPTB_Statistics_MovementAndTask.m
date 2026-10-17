function [Statistics] = GSPTB_Statistics_MovementAndTask(Measure)

    assert(isnumeric(Measure) && ismatrix(Measure), ...
    'Measure must be a subjects-by-tasks numeric matrix.');

    assert(size(Measure,2) > 1, ...
    'Measure must contain at least two tasks.');

    % Number of subjects
    S = size(Measure,1);

    % Number of tasks
    n_tasks = size(Measure,2);

    % Conducts repeated-measures modeling
    VariableNames = compose("Task%d",1:n_tasks);
    Data = array2table(Measure,'VariableNames',VariableNames);
    
    WithinDesign = table(categorical((1:n_tasks)'), ...
        'VariableNames',{'Task'});
    
    rm = fitrm(Data,sprintf('%s-%s ~ 1',VariableNames(1),VariableNames(end)), ...
        'WithinDesign',WithinDesign);
    
    % Subsequent repeated-measures ANOVA
    aov_rm = ranova(rm,'WithinModel','Task');

    % Construct the group labels
    Group_Labels = 1:n_tasks;
    Group_Labels = repmat(Group_Labels,S,1);
    
    % Conducts a simple ANOVA to assess the measure of interest
    % aov = anova(Group_Labels(:),Measure(:));

    % Initializes the statistical arrays
    p_values = NaN(n_tasks,n_tasks);
    t_stats = NaN(n_tasks,n_tasks);
    Cohen_d = NaN(n_tasks,n_tasks);
    Cohen_dz = NaN(n_tasks,n_tasks);

    % Statistics used for subsequent Cohen's d quantifications
    Mu = mean(Measure);
    Sigma2 = var(Measure);
    
    % Computes individual contrasts
    for t1 = 1:n_tasks
        for t2 = 1:n_tasks
            if t2 ~= t1

                % p-values
                [~,p_values(t1,t2),~,stats] = ttest(Measure(:,t1),Measure(:,t2));

                % t-statistics
                t_stats(t1,t2) = stats.tstat;

                % Cohen's d_z
                Difference = Measure(:,t1) - Measure(:,t2);
                Cohen_dz(t1,t2) = mean(Difference,'omitnan') / std(Difference,'omitnan');

                % Cohen's d
                Cohen_d(t1,t2) = (Mu(t1)-Mu(t2))/sqrt((Sigma2(t1)+Sigma2(t2))/2);
            end
        end
    end

    % Now performs a linear mixed effects model, where subject is modeled
    % as a random intercept
    Subject_vec = repelem((1:S)',n_tasks);
    Task_vec = repmat((1:n_tasks)',S,1);
    Measure_vec = reshape(Measure.',[],1);

    DataTab = table(categorical(Subject_vec),categorical(Task_vec),Measure_vec,...
    'VariableNames',{'Subject','Task','Measure'});

    lme = fitlme(DataTab,'Measure ~ Task + (1|Subject)',...
    'FitMethod','REML');

    % We extract the residual and subject variances and compute the ICC
    [Psi,MSE] = covarianceParameters(lme);

    Variance_Subject = Psi{1}(1,1);
    Variance_Residual = MSE;

    ICC = Variance_Subject/(Variance_Subject + Variance_Residual);

    % Fills in the structure of statistics to return
    % Statistics.ANOVA = aov;
    Statistics.RepeatedMeasuresANOVA = aov_rm;
    Statistics.p_values = p_values;
    Statistics.t_statistics = t_stats;
    Statistics.Cohen_d = Cohen_d;
    Statistics.Cohen_dz = Cohen_dz;
    Statistics.ICC = ICC;
end