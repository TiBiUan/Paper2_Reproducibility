function [X_OI] = GSPTB_LoadTaskData(FuncPath,Task_name,idx_cortical,is_cortical)

    Current = pwd;

    % We go and load the functional data at hand
    cd(FuncPath);

    name = ['X_',Task_name,'.mat'];
    Loaded = load(name);
    VariableName = ['X_',Task_name];
    X_OI = Loaded.(VariableName);

    S = length(X_OI);

    % We must remove subcortical regions if doing cortex-only analyses
    if is_cortical
        for s = 1:S
            for ses = 1:2
                X_OI{s}{ses} = X_OI{s}{ses}(:,idx_cortical);
            end
        end
    end
 
    cd(Current);
end