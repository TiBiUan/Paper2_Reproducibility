function [G] = GSPTB_ComputePopulationConsensusMask(A,idx,dist,hemiid,nbins)

    if isempty(idx)
        error('Vector idx should not be empty!');
    end
    
    % Reformats selected SC matrices into R x R x N format
    A_Selected = cat(3,A{idx});
    
    % Gets the binary distance-dependent consensus mask
    G = fcn_group_bins(A_Selected,dist,hemiid,nbins);
end