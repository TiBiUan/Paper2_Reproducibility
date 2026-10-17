function [Results] = GSPTB_LME_TaskContrasts(Models,TaskVariable,Alpha)

    % We always use Satterthwaite degree of freedom as we work with mixed
    % models with a random subject effect
    DFMethod = 'satterthwaite';

    % Samples the information about data points in terms of their task tag
    TaskData = Models{1}.Variables.(TaskVariable);

    % All task levels
    TaskLevels = categories(TaskData);

    % Number of tasks
    nTasks = numel(TaskLevels);

    % Number of harmonics (one per model)
    nHarmonics = numel(Models);
    
    % Initialization for outputs
    MatrixSize = [nHarmonics,nTasks,nTasks];

    Estimate = NaN(MatrixSize);
    SE = NaN(MatrixSize);
    tStat = NaN(MatrixSize);
    FStat = NaN(MatrixSize);
    DF1 = NaN(MatrixSize);
    DF2 = NaN(MatrixSize);
    pValue = NaN(MatrixSize);
    CI_Lower = NaN(MatrixSize);
    CI_Upper = NaN(MatrixSize);

    % List of unique comparisons
    PairIndices = nchoosek(1:nTasks,2);
    nPairs = size(PairIndices,1);

    ContrastLabels = strings(nPairs,1);

    for k = 1:nPairs
        ContrastLabels(k) = string(TaskLevels{PairIndices(k,1)}) + " - " + ...
            string(TaskLevels{PairIndices(k,2)});
    end

    % Same process for each harmonic...
    for c = 1:nHarmonics

        disp(['Harmonic ',num2str(c),'...']);
        
        % Samples the current model
        lme = Models{c};

        CurrentTask = lme.Variables.(TaskVariable);
       
        % Fixed effect quantities

        % Design matrix
        X = designMatrix(lme,'Fixed');

        % Fixed effect coefficients
        Beta = fixedEffects(lme);

        % Covariance across coefficients
        % This tells us the relationship between the different coefficients
        % around the returned optimum
        VBeta = lme.CoefficientCovariance;

        % Name of all considered coefficients
        CoefficientNames = string(lme.CoefficientNames);

        % Number of coefficients
        nCoefficients = numel(Beta);

        % Specific 
        XTask = NaN(nTasks,nCoefficients);

        for t = 1:nTasks

            % Only flags data for task t
            idx = (CurrentTask == TaskLevels{t});

            % Finds the index of the first time the task at hand happens
            firstRow = find(idx,1,'first');
            XTask(t,:) = X(firstRow,:);
        end

        % In the intended model, the only continuous covariate is Motion.
        % Setting its fixed-effect design column to zero evaluates all task
        % means at the centered Motion value.
        %
        % More importantly, every task row receives the same value, so the
        % Motion term cancels exactly in each task contrast.
        
        % Sets the motion column of XTask to 0. This is because we want to
        % examine task contrasts independently from motion
        MotionColumn = find(strcmpi(CoefficientNames,'Motion'));

        if numel(MotionColumn) ~= 1
            error(['Expected exactly one fixed-effect coefficient named Motion, ', ...
                   'but found %d. Coefficients are: %s'], ...
                   numel(MotionColumn), ...
                   strjoin(CoefficientNames,', '));
        end
       
        % We consider all individual task contrasts...
        for k = 1:nPairs

            % The two tasks that must be compared
            t1 = PairIndices(k,1);
            t2 = PairIndices(k,2);

            % Contrast orientation:
            %
            %     Task t1 minus Task t2
            H = XTask(t1,:) - XTask(t2,:);

            % Adjusted task difference.
            CurrentEstimate = H*Beta;

            % We compute the variance of our contrast estimate (how
            % much it may vary if we repeated the process with another
            % dataset)
            CurrentVariance = H*VBeta*H';

            Tolerance = 1e-10 * max(1,norm(VBeta,'fro') * norm(H)^2);

            if CurrentVariance < -Tolerance
                error('Substantially negative contrast variance: %.12g', ...
                      CurrentVariance);
            elseif CurrentVariance < 0
                CurrentVariance = 0;
            end

            CurrentSE = sqrt(CurrentVariance);

            % Inferential test.
            [CurrentP,CurrentF,CurrentDF1,CurrentDF2] = ...
                coefTest(lme,H,0,'DFMethod',DFMethod);

            % A one-dimensional contrast satisfies F = t^2.
            CurrentT = sign(CurrentEstimate)*sqrt(CurrentF);

            % Confidence interval.
            CriticalValue = tinv(1-Alpha/2,CurrentDF2);
            
            Lower = CurrentEstimate - CriticalValue*CurrentSE;
            Upper = CurrentEstimate + CriticalValue*CurrentSE;

            %% Store the requested orientation: t1 - t2

            Estimate(c,t1,t2) = CurrentEstimate;
            SE(c,t1,t2)       = CurrentSE;
            tStat(c,t1,t2)    = CurrentT;
            FStat(c,t1,t2)    = CurrentF;
            DF1(c,t1,t2)      = CurrentDF1;
            DF2(c,t1,t2)      = CurrentDF2;
            pValue(c,t1,t2)   = CurrentP;
            CI_Lower(c,t1,t2) = Lower;
            CI_Upper(c,t1,t2) = Upper;

            %% Store the reverse orientation: t2 - t1

            Estimate(c,t2,t1) = -CurrentEstimate;
            SE(c,t2,t1)       = CurrentSE;
            tStat(c,t2,t1)    = -CurrentT;
            FStat(c,t2,t1)    = CurrentF;
            DF1(c,t2,t1)      = CurrentDF1;
            DF2(c,t2,t1)      = CurrentDF2;
            pValue(c,t2,t1)   = CurrentP;

            % Reversing the estimate also reverses and swaps the CI bounds.
            CI_Lower(c,t2,t1) = -Upper;
            CI_Upper(c,t2,t1) = -Lower;
        end
    end

    %% Returning the results
    Results = struct;

    Results.Estimate = Estimate;
    Results.SE = SE;
    Results.tStat = tStat;
    Results.FStat = FStat;
    Results.DF1 = DF1;
    Results.DF2 = DF2;
    Results.pValue = pValue;
    Results.CI_Lower = CI_Lower;
    Results.CI_Upper = CI_Upper;

    Results.TaskVariable = TaskVariable;
    Results.TaskLevels = TaskLevels;
    Results.PairIndices = PairIndices;
    Results.ContrastLabels = ContrastLabels;

    Results.DFMethod = DFMethod;
    Results.Alpha = Alpha;
end