function [Results] = GSPTB_RunHarmonicWisePredictionCurve(U,Lambda,DegMat,...
    PredType,b,PredOpts)

    % Number of harmonics
    R = size(U{1},1);

    % Initialization
    R2_SC = NaN(R,PredOpts.NumOuterSplits);
    R2_SC_cova = NaN(R,PredOpts.NumOuterSplits);
    
    % We consider all harmonics sequentially...
    for c = 2:R

        disp(['Currently considering case number ',num2str(c),'...']);

        % Specifies the range of harmonics to use for prediction
        switch PredType
            case 'cumlow'
                SelRange = 1:c;
            case 'cumhigh'
                SelRange = (R-c+1):R;
            case 'indiv'
                SelRange = c;
        end
    
        Ker_SC = GSPTB_ComputeHarmonicRangeKernel(U,Lambda,DegMat,...
            PredOpts.NormalizeHarmonicAdjacency,SelRange);
        
        % Performs the prediction itself
        PredRes = GSPTB_KernelRidgeRegression(Ker_SC,...
            b,PredOpts);

        R2_SC(c,:) = PredRes.R2;
        R2_SC_cova(c,:) = PredRes.R2_Cova;
    end

    Results.NoCovariate = R2_SC;
    Results.Covariate = R2_SC_cova;
end