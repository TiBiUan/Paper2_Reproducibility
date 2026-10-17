function [PredictionResults] = GSPTB_RunSingleHarmonicWisePrediction(U,Lambda,DegMat,...
    Harm_Range,b,PredOpts)

    % Computes the structural kernel for the harmonic range at hand
    Ker_SC = GSPTB_ComputeHarmonicRangeKernel(U,Lambda,DegMat,...
            PredOpts.NormalizeHarmonicAdjacency,Harm_Range);
    
    % Performs the prediction itself
    PredictionResults = GSPTB_KernelRidgeRegression(Ker_SC,b,PredOpts);
end