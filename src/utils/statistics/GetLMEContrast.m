function Out = GetLMEContrast(lme,C)

    beta  = fixedEffects(lme);
    Vbeta = lme.CoefficientCovariance;

    Estimate = C*beta;
    SE = sqrt(C*Vbeta*C');

    [p,~,~,DF2] = coefTest(lme,C);

    tcrit = tinv(0.975,DF2);

    CI_Lower = Estimate-tcrit*SE;
    CI_Upper = Estimate+tcrit*SE;

    Out = table(Estimate,SE,CI_Lower,CI_Upper,p);
end