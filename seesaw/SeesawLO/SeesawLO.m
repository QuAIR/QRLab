function [a, choi_A, choi_B] = SeesawLO(JN_target, dA, dB, varargin)
    
    % .. math::
    %
    %     \gamma_{LO}(\mathcal{N}_{AB \rightarrow A'B'}) = \min \sum_j |\alpha_j| \
    %     s.t. \
    %     \mathcal{M}^{(A/B)}_j \in \operatorname{CPTN}(A/B \rightarrow A'/B') \
    %     \mathcal{N}_{AB\rightarrow A'B'} = \sum_j \alpha_j \mathcal{M}^{(A)}_j \otimes \mathcal{M}^{(B)}_j     
    %
    % Args:
    %   JN_target (numeric): The Choi matrix of the bipartite channel.
    %   dA (numeric): The dimension of subsystem A.
    %   dB (numeric): The dimension of subsystem B.
    %   varargin (numeric): The hyperparameters for the optimization algorithm.
    %
    % Returns:
    %   numeric: The optimized decomposition of the target channel.

    % optional parameters
    p = inputParser;
    addRequired(p, 'JN_target');
    addRequired(p, 'dA');
    addRequired(p, 'dB');
    addOptional(p, 'num_pos', 20);
    addOptional(p, 'num_neg', 20);
    addOptional(p, 'num_itr', 400);
    addOptional(p, 'error_criteria', 1e-7);
    addOptional(p, 'cost_convergence_criteria', 5e-7);
    % addOptional(p, 'convergence_criteria', 50);
    addOptional(p, 'c', 0.9999);
    
    parse(p, JN_target, dA, dB, varargin{:});
    
    % extract parameters
    JN_target = p.Results.JN_target;
    dA = p.Results.dA;
    dB = p.Results.dB;
    num_pos = p.Results.num_pos;
    num_neg = p.Results.num_neg;
    num_term = num_pos + num_neg;
    num_itr = p.Results.num_itr;
    error_criteria = p.Results.error_criteria;
    cost_convergence_criteria = p.Results.cost_convergence_criteria;
    % convergence_criteria = p.Results.convergence_criteria;
    c = p.Results.c;

    % distortion
    distortion = PermuteSystems(kron(RandomSuperoperator(dA), RandomSuperoperator(dB)), [1,3,2,4], [dA, dA, dB, dB]);
    distort_target = (1-c)*JN_target + c*distortion;

    %% Seesaw decomposition
    fprintf('\n Seesaw decomposition start! \n');
    %TODO: Generalize to SEP, LOCC, etc.
    % initialization A
    choi_A = zeros(dA^2, dA^2, num_term);
    for j=1:num_term
        choi_A(:,:,j) = RandomSuperoperator(dA);
    end

    err_list = {};
    for idx=1:num_itr

        % optimize B
        [a, choi_B, err] = SeesawDecoupling(JN_target, choi_A, num_pos, 2);
        err_list{end+1} = err;

        % optimize A
        [a, choi_A, err] = SeesawDecoupling(JN_target, choi_B, num_pos, 1);
        err_list{end+1} = err;

        if err < error_criteria
            break;
        elseif numel(err_list) > 100 && abs(err_list{end} - err_list{end-1}) < error_criteria
            break;
        else
            fprintf('The minimum decoupling error does not meet the convergence criteria: received %d \n', err);
            return;
        end
    end
    a = a; choi_A = choi_A; choi_B = choi_B; err = err;
    fprintf('\n The minimum decoupling error reaches %d\n', err);

    fprintf('Seesaw cost minimization start! \n');
    %% Seesaw minimization
    %TODO: Generalize to SEP, LOCC, etc.

    cost_list = {};
    for idx=1:num_itr
        
        % optimize B
        [a, choi_B] = SeesawMinimizing(JN_target, choi_A, error_criteria, num_pos, 2);
        cost_list{end+1} = sum(a);

        % optimize A
        [a, choi_A] = SeesawMinimizing(JN_target, choi_B, error_criteria, num_pos, 1);
        cost_list{end+1} = sum(a);

        if cost_list{end} < 1.001 || abs(cost_list{end} - cost_list{end-1}) < cost_convergence_criteria
            break;
        end
    end
    fprintf('The minimum cost reaches %d \n', cost_list{end});

    % collecting results
    a = a;
    choi_A = choi_A;
    choi_B = choi_B;
end