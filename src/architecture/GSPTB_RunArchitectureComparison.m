function [Kernel,ArchitectureResults] = GSPTB_RunArchitectureComparison(U_R,Power,b,Cfg,PredOpts)

    % Number of subjects
    S = length(U_R);

    KernelTypes = {'Structure','Function','Hybrid'};
    PairingTypes = {'Matched','All'};

    % Generates the permutations; we will use the same across all cases to
    % minimize sources of variability
    id_perm = NaN(Cfg.Statistics.Prediction.NumNulls,S);

    for n = 1:Cfg.Statistics.Prediction.NumNulls
        id_perm(n,:) = randperm(S);
    end
    
    for i = 1:length(KernelTypes)
        for j = 1:length(PairingTypes)

            % Computes the kernel to study
            Kernel.(KernelTypes{i}).(PairingTypes{j}) = GSPTB_CreateKernel(U_R,Power(Cfg.Tasks.MainFigureIndices,:,:),KernelTypes{i},PairingTypes{j});

            % Performs the prediction itself
            PredResults = GSPTB_KernelRidgeRegression(Kernel.(KernelTypes{i}).(PairingTypes{j}),b,PredOpts);

            ArchitectureResults.(KernelTypes{i}).(PairingTypes{j}) = PredResults;
            ArchitectureResults.(KernelTypes{i}).(PairingTypes{j}).NullShuffles = id_perm;

            % Now performs the null prediction
            % for n = 1:Cfg.Statistics.Prediction.NumNulls
            % 
            %     PredResultsNull = GSPTB_KernelRidgeRegression(Kernel.(KernelTypes{i}).(PairingTypes{j}),b(id_perm(n,:)),PredOpts);
            %     ArchitectureResults.(KernelTypes{i}).(PairingTypes{j}).NullData{n} = PredResultsNull;
            % end
        end 
    end
end