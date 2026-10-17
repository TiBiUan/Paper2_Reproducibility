function [MovementStatistics] = GSPTB_RunMovementStatistics(Data)

    % Establishes movement statistics across tasks for all motion
    % descriptions we have in the data
    MovementStatistics.MeanFD = GSPTB_Statistics_MovementAndTask(Data.MeanFD);
    MovementStatistics.MeanFD_NonScrubbedOnly = GSPTB_Statistics_MovementAndTask(Data.MeanFD_NonScrubbedOnly);
    MovementStatistics.PercentageScrubbed = GSPTB_Statistics_MovementAndTask(Data.PercentageScrubbed);
end