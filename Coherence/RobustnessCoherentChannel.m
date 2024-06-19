function C_R = RobustnessCoherentChannel(channel, varargin)
    % Provide the robustness of a coherent channel.
    %
    % .. math::
    %    1 + C_R(\mathcal{N}) = \min_{\mathcal{M}}\{\lambda | \mathcal{N} \le \lambda \mathcal{M}\},
    %
    % where :math:`\mathcal{M}` is a maximally incoherent operation (MIO).
    %
    % Args:
    %     channel (matrix): The Choi matrix of the input channel :math:`\mathcal{N}`.
    %     varargin (list): Dimension of the given channel, default to
    %      [d_i,d_o], with d_i=d_o.
    %
    % Returns:
    %     C_R (numeric): The robustness of the input channel.
    %
    % Note:
    %     Díaz, M. G., Fang, K., Wang, X., Rosati, M., Skotiniotis, M., Calsamiglia, J., & Winter, A. (2018). 
    %     Using and reusing coherence to realize quantum processes. 
    %     Quantum, 2, 100.


% set optional argument defaults
dim = opt_args({[sqrt(length(channel)), sqrt(length(channel))]}, varargin{:});
d_i = dim(1);
d_o = dim(2);
JN = channel;

cvx_begin sdp quiet

    variable JM(d_i*d_o, d_i*d_o) hermitian
    variable s

    minimize s
    subject to
        JM >= 0;
        PartialTrace(JM,2, [d_i,d_o]) == (1+s)*eye(d_i)
        s >= 0;
        JN <= JM;
        % MIO constraints
        for i = 1:d_i
        PartialTrace(JM*kron(KetBra(d_i,i,i), eye(d_o)),1,[d_i, d_o])...
            == diag(diag(PartialTrace(JM*kron(KetBra(d_i,i,i), eye(d_o)),1,[d_i, d_o])));
        end
cvx_end
C_R = s;

end