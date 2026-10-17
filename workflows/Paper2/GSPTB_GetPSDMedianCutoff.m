function [id_low,id_high] = GSPTB_GetPSDMedianCutoff(X_TC,FD,H,id,session_id,T_FD)

    % Builds the data for all subjects to consider
    for s = 1:numel(id)

        % Time courses of resting-state activity
        tmp = X_TC{id(s)};
        X_OI{s} = tmp{session_id};

        % Time course of framewise displacement
        FD_OI{s} = FD{id(s)};

        % Gets the spectral representation
        Xhat_OI{s} = GSPTB_GFT(H.U,X_OI{s}(FD_OI{s}<T_FD,:)');

        % PSD
        PSD(:,s) = GSPTB_ComputePSD(Xhat_OI{s});
    end

    % Mean PSD across the population
    MU_PSD = mean(PSD,2);

    % Finds the median cutoff
    CPSD = cumsum(MU_PSD);

    id_cutoff = find(CPSD > 0.5*max(CPSD));
    id_cutoff = id_cutoff(1);

    id_low = 1:id_cutoff-1;
    id_high = id_cutoff:length(CPSD);
end