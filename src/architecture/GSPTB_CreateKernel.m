function [Ker] = GSPTB_CreateKernel(U,FuncMeasure,KernelType,PairingType)

    S = length(U);
    R = size(U{1},1);

    % initializes the kernel
    Ker = NaN(S,S);

    % We compute our kernel by looking at all subject pairs...
    for s1 = 1:S

        s1
        
        for s2 = s1:S
            switch KernelType
                case 'Structure'

                    % Structural similarity weights
                    W = abs(corr(U{s1},U{s2}));

                case 'Function'

                    % Functional similarity weights
                    W = corr(squeeze(FuncMeasure(:,s1,:)),squeeze(FuncMeasure(:,s2,:)));
                
                case 'Hybrid'
                    
                    % Both subtypes of weights combined
                    W_S = abs(corr(U{s1},U{s2}));
                    W_F = corr(squeeze(FuncMeasure(:,s1,:)),squeeze(FuncMeasure(:,s2,:)));
                    W = W_S.*W_F;
                otherwise
                    error('Unknown kernel type!');
            end

            switch PairingType
                case 'Matched'
                    Ker(s1,s2) = mean(diag(W));
                    
                case 'All'
                    Ker(s1,s2) = mean(W,'all'); 
                otherwise
                    error('Unknown pairing type!');
            end

            % Fills the lower triangle
            Ker(s2,s1) = Ker(s1,s2);
        end
    end
end

