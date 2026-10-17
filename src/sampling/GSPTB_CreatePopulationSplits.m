% Splits, aware of family relationships, between left-out and retained data
% Then, on the retained data, sample a subset of subjects; have increasing
% sample sizes with nested subject selection
function [PopSplits] = GSPTB_CreatePopulationSplits(Opts)

    % Making sure the values of sample sizes are monotonically increasing
    if any(diff(Opts.SampleSizes) <= 0)
        error('Opts.SampleSizes must contain strictly increasing sample sizes!');
    end

    % Making sure sample sizes are positive integers
    if any(Opts.SampleSizes <= 0) || any(mod(Opts.SampleSizes,1) ~= 0)
        error('Opts.SampleSizes must contain only positive integer values!');
    end

    % Number of subjects
    S = size(Opts.FamilyKernel,1);
    
    % How many subjects in each side of the outer split
    n_leftout = floor(Opts.LeftOutFraction*S);
    n_tt = S - n_leftout;

    % Splitting for each outer fold
    for no = 1:Opts.NumOuterSplits
    
        % Splitting (outer loop) keeping track of family structure
        [idx_tt{no},idx_leftout{no}] = GSPTB_CreateSplit(n_tt,...
            n_leftout,Opts.FamilyKernel);

        % Making sure all requested sample sizes are compatible with the
        % actual retained subject pool for this outer split
        if any(Opts.SampleSizes > numel(idx_tt{no}))
            error(['Opts.SampleSizes cannot contain values larger than the ' ...
                'number of retained subjects in the current outer split!']);
        end
        
        % The inner process of sampling is performed NumMonteCarlo times
        for nm = 1:Opts.NumMonteCarlo

            % Current shuffling of indices to consider
            perm = randperm(numel(idx_tt{no}));

            % Just for clearer following code lines
            order = idx_tt{no}(perm);

            % Full ordering for MC fold nm; note that we want the indices
            % relative to the full cohort, so we shuffle idx_tt{no}
            PopSplits.Outer(no).MonteCarlo(nm).Order = order;

            % We sample increasing numbers of subjects from the same ordering,
            % thereby generating nested samples
            for nn = 1:numel(Opts.SampleSizes)
        
                % Current number of samples
                N = Opts.SampleSizes(nn);
            
                % Gets the actual samples for the current outer loop index
                % (no), the current Monte Carlo sampling index (nm), and
                % the current number of samples (nn)
                PopSplits.Outer(no).MonteCarlo(nm).SampleA{nn} = order(1:N);

                % If we consider a small enough current sample size that
                % two disjoint sets can be built, we create the second set
                % (used for robustness analyses)
                if 2*N <= numel(order)
                    PopSplits.Outer(no).MonteCarlo(nm).SampleB{nn} = order(end-N+1:end);

                % Else, we just fill with NaN instead
                else
                    PopSplits.Outer(no).MonteCarlo(nm).SampleB{nn} = [];
                end
            end
        end

        % Indices of training (undergoing MC sampling) and test (left-out,
        % to be used for representativeness) subsets
        PopSplits.Outer(no).Train = idx_tt{no};
        PopSplits.Outer(no).Test = idx_leftout{no};

        % Associated exact numbers of samples
        PopSplits.Outer(no).NumTrain = numel(idx_tt{no});
        PopSplits.Outer(no).NumTest = numel(idx_leftout{no});
    end

    % Saving sample sizes as well, and other important inputs
    PopSplits.SampleSizes = Opts.SampleSizes;
    PopSplits.NumOuterSplits = Opts.NumOuterSplits;
    PopSplits.NumMonteCarlo = Opts.NumMonteCarlo;
    PopSplits.LeftOutFraction = Opts.LeftOutFraction;
end