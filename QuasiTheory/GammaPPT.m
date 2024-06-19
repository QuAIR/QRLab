function val = GammaPPT(JN, dim)

    % .. math::
    %
    %     \gamma_{\operatorname{PPT}}(\mathcal{N}) = \min\{p_+ + p_-:\mathcal{N} = p_+\mathcal{M_+} - p_-\mathcal{M_-}\},
    %
    % Args:
    %   JN (numeric): The Choi matrix of the bipartite channel.
    %   dim (numeric): The array storing input and output dimensions.
    %
    % Returns:
    %   numeric: The PPT-assisted sampling overhead of bipartite channel.
    %
    % Raises:
    %   error: If either input/output dimension does not match, an error is raised.
    %
    % Note:
    %   Jing, M., Zhu, C., & Wang, X. (2024). 
    %   Circuit Knitting Faces Exponential Sampling Overhead Scaling Bounded by Entanglement Cost. 
    %   arXiv preprint arXiv:2404.03619.

    dA = dim(1);
    dB = dim(2);
    dAB = dA*dB;

    cvx_begin sdp quiet
    variable JM(dAB^2, dAB^2) hermitian
    variable c nonnegative
        
    cost = 2*c - 1;
    minimize cost
    subject to
        JM >= JN; PartialTranspose(JM, [2, 4], [dA, dB, dA, dB]) >= 0;
        PartialTrace(JM, 2, [dAB dAB]) == c * eye(dAB);
        PartialTranspose(JM, [2, 4], [dA, dB, dA, dB]) >= PartialTranspose(JN, [2, 4], [dA, dB, dA, dB]);
    cvx_end 
    
    val = cost;
end