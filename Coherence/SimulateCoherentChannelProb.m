function probability = SimulateCoherentChannelProb(target_channel, resource_state, free_op, error_tolerance, varargin)
    % Find the maximally success probability :math:`p` to simulate the target channel
    % simulation with free operations :math:`\mathcal{M}` and resource state
    % :math:`\omega` up to error tolerance :math:`\epsilon`.
    %
    % .. math::
    %    \max_{\mathcal{M}} \{p| \mathcal{M}(\sigma\otimes\cdot)=p\mathcal{L}(\cdot) 
    %    + (1-p)\mathcal{R}(\cdot), \|\mathcal{L} - \mathcal{N}\|\le \epsilon\}
    %
    % where :math:`\mathcal{R}` is a rubbish channel.
    %
    % Args:
    %     target_channel (matrix): The Choi matrix of the input channel :math:`\mathcal{N}`.
    %     reousrce_state (matrix): The density matrix of the given resource
    %      state :math:`\omega`.
    %     free_op (numeric): The choice of free operation
    %      :math:`\mathcal{M}`, we can choose `0` (MIO) or `1` (DIO).
    %     error_tolerance (numeric): The maximal simulation error allowed.
    %     varargin (list): Dimension of the given channel, default to
    %      `[d_i,d_o]`, with `d_i=d_o`.
    %
    % Returns:
    %     probability (numeric): The maximally success probability of channel
    %     simulation.
    %
    % Note:
    %     Zhao, B., Ito, K., & Fujii, K. (2024). 
    %     Probabilistic channel simulation using coherence. 
    %     arXiv preprint arXiv:2404.06775.

% set optional argument defaults
dim = opt_args({[sqrt(length(target_channel)), sqrt(length(target_channel))]}, varargin{:});
dim_i = dim(1);
dim_o = dim(2);
dim_r = length(resource_state);
JN = target_channel;
eps = error_tolerance;
%% SDP
cvx_begin sdp quiet
    variable JM(dim_r*dim_i*dim_o,dim_r*dim_i*dim_o) hermitian
    variable Z(dim_i*dim_o,dim_i*dim_o) hermitian
    variable t

    JE = PartialTrace(JM*kron(resource_state.', eye(dim_i*dim_o)), 1, [dim_r,dim_i*dim_o]);
    minimize t
    subject to
        JM >= 0;
        PartialTrace(JM, 3, [dim_r,dim_i,dim_o]) <= t * eye(dim_r*dim_i);
        PartialTrace(JE, 2, [dim_i,dim_o]) == eye(dim_i);
        Z >= 0;
        Z >= JN-JE;
        eps*eye(dim_i) >= PartialTrace(Z,2,[dim_i, dim_o]);
        %%%%% MIO %%%%%%%%%
        for i = 1:dim_r*dim_i
            PartialTrace(JM*kron(KetBra(dim_r*dim_i,i,i), eye(dim_o)),1,[dim_r*dim_i, dim_o])...
                == diag(diag(PartialTrace(JM*kron(KetBra(dim_r*dim_i,i,i), eye(dim_o)),1,[dim_r*dim_i, dim_o])));
        end

        if free_op == 1
            %%%%% DIO %%%%%%%%%
            for i = 1:dim_r*dim_i
                for j = 1:dim_r*dim_i
                    if i ~= j
                        diag(diag(PartialTrace(JM*kron(KetBra(dim_r*dim_i, i, j).', eye(dim_o)), 1, [dim_r*dim_i,dim_o]))) == zeros(dim_o,dim_o);
                    end
                end
            end
        end
cvx_end
probability = 1/t;
