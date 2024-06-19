function distance = SimulateCoherentChannel(target_channel, resource_state, free_op, varargin)
    % Find the optimal free operation :math:`\mathcal{M}` to minimize
    % diamond norm between the simulated channel
    % :math:`\mathcal{M}(\omega\otimes\cdot)` and the target channel
    % :math:`\mathcal{N}` 
    %
    % .. math::
    %    \| \mathcal{M}(\omega \otimes \cdot) - \mathcal{N}\|_\diamond,
    %
    % where :math:`\omega` is the given resource state.
    %
    % Args:
    %     target_channel (matrix): The Choi matrix of the input channel :math:`\mathcal{N}`.
    %     reousrce_state (matrix): The density matrix of the given resource
    %      state :math:`\omega`.
    %     free_op (numeric): The choice of free operation
    %      :math:`\mathcal{M}`, we can choose 0 (MIO) or 1 (DIO).
    %     varargin (list): Dimension of the given channel, default to
    %      [d_i,d_o], with d_i=d_o.
    %
    % Returns:
    %     distance (numeric): The diamond norm between the simulated channel and the target channel. 
    %
    % Note:
    %     Díaz, M. G., Fang, K., Wang, X., Rosati, M., Skotiniotis, M., Calsamiglia, J., & Winter, A. (2018). 
    %     Using and reusing coherence to realize quantum processes. 
    %     Quantum, 2, 100.



% set optional argument defaults
dim = opt_args({[sqrt(length(target_channel)), sqrt(length(target_channel))]}, varargin{:});
dim_i = dim(1);
dim_o = dim(2);
dim_r = length(resource_state);
JN = target_channel;

%% SDP
cvx_begin sdp quiet
    variable JM(dim_r*dim_i*dim_o,dim_r*dim_i*dim_o) hermitian
    variable Z(dim_i*dim_o,dim_i*dim_o) hermitian
    variable l

    JE = PartialTrace(JM*kron(resource_state, eye(dim_i*dim_o)), 1, [dim_r,dim_i,dim_o]);
    minimize l
    subject to
        JM >= 0;
        PartialTrace(JM, 3, [dim_r,dim_i,dim_o]) == eye(dim_r*dim_i);

        Z >= 0;
        Z >= JN-JE;
        l*eye(dim_i) >= PartialTrace(Z,2, [dim_i, dim_o]);

        PartialTrace(JE,2,[dim_i,dim_o]) == eye(dim_i);
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
distance = l;