% Computes similarity between two harmonic bases
function [Sim] = GSPTB_ComputeHarmonicSimilarity(U1,U2,match_type)

    if ~isequal(size(U1),size(U2))
        error('U1 and U2 must have the same dimensions!');
    end

    switch match_type
        case 'index'
            Sim = diag(abs(corr(U1,U2)));
        case 'match'
            % Aligns U2 to U1, then computes similarity afterwards
            U2_aligned = GSPTB_RealignHarmonics_Hungarian(U2,U1,zeros(size(U1,1),1));
            Sim = diag(abs(corr(U1,U2_aligned)));
        otherwise
            error('Subtype of matching not found!');
    end
end