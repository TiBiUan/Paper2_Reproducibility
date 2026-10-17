function [R2_Null,id_perm] = GSPTB_RunNullHarmonicWisePrediction(U,...
    Lambda,DegMat,Harm_Range,b,n_null_behavpred,is_cova,PredOpts)

    % Number of subjects
    S = length(U);

    % Generates the permutations for the behavioral score
    id_perm = NaN(n_null_behavpred,S);
    
    % Initializes the R2 arrays
    R2_Null = NaN(n_null_behavpred,PredOpts.NumOuterSplits);
    
    % Constructs the kernel (identical for all nulls)
    Ker_SC = GSPTB_ComputeHarmonicRangeKernel(U,Lambda,DegMat,...
                PredOpts.NormalizeHarmonicAdjacency,Harm_Range);
    
    % Performs the null permutations themselves
    for n = 1:n_null_behavpred
        
        % We make sure to save the permutations
        id_perm(n,:) = randperm(S);
    
        PredRes = GSPTB_KernelRidgeRegression(Ker_SC,b(id_perm(n,:)),PredOpts);

        if ~is_cova
            R2_Null(n,:) = PredRes.R2;
        else
            R2_Null(n,:) = PredRes.R2_Cova;
        end
    end
end