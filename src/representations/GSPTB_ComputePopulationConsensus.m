function [A_Pop] = GSPTB_ComputePopulationConsensus(A,idx,dist,hemiid,nbins)

    % SC matrices to consider
    A_touse = A(idx);

    % Reformatted to be compatible
    for s = 1:length(A_touse)
        A_OI(:,:,s) = A_touse{s};
    end
    
    % Gets the binary consensus
    Consensus_Binary = fcn_group_bins(A,dist,hemiid,nbins);
    Weights = squeeze(mean(A_OI,3));

    % Gets the output
    A_Pop = Weights .* Consensus_Binary;
end