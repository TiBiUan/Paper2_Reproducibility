function [BehavioralPrediction] = GSPTB_RunFamilyWisePrediction(U_R,Lambda_R,DegMat,idx_fam,b,PredOpts)

    % Number of families
    F = max(idx_fam);
    
    % Now, we want to perform family-wise prediction
    R2_SingleFam = NaN(max(idx_fam),PredOpts.NumOuterSplits);
    R2_CumLowFam = NaN(max(idx_fam),PredOpts.NumOuterSplits);

    % We consider all harmonic families sequentially
    for fam = 1:F

        PredRes = GSPTB_RunSingleHarmonicWisePrediction(U_R,Lambda_R,DegMat,find(idx_fam==fam),b,PredOpts);
        R2_SingleFam(fam,:) = PredRes.R2;
    
        PredRes = GSPTB_RunSingleHarmonicWisePrediction(U_R,Lambda_R,DegMat,find(idx_fam<=fam),b,PredOpts);
        R2_CumLowFam(fam,:) = PredRes.R2;
    end

    BehavioralPrediction.Single = R2_SingleFam;
    BehavioralPrediction.CumLow = R2_CumLowFam;
    
    Delta_R2_CumLowFam = diff(R2_CumLowFam);
    Delta_R2_CumLowFam = [R2_SingleFam(1,:);Delta_R2_CumLowFam];

    BehavioralPrediction.Incremental = Delta_R2_CumLowFam;
end