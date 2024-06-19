function val = MaxLogNeg(JN, dim)
    % .. math::
    %
    %     LN_{\max}(\mathcal{N}) = \log \inf\{\max\{\|P_{AB}\|_{\infty},
    %     \|P^{T_B}_{AB}\|_{\infty}\}: -P^{T_{BB'}}_{ABA'B'} \leq (J^{\mathcal{N}}_{ABA'B'})^{T_{BB'}} \leq P^{T_{BB'}}_{ABA'B'}\},
    %
    % Args:
    %     JN (numeric): The Choi matrix of the bipartite channel.
    %     dim (numeric): The array storing input and output dimensions.
    %
    % Returns:
    %     numeric: The generalized :math:`\kappa`-entanglement 
    %     (or max-logarithmic negativity) of bipartite channel.
    %
    % Raises:
    %     error: If either input/output dimension does not match, an error is raised.
    %
    % Note:
    %     Wang, X., & Wilde, M. M. (2023). 
    %     Exact entanglement cost of quantum states and channels under positive-partial-transpose-preserving operations. 
    %     Physical Review A, 107(1), 012429.


    dA = dim(1);
    dB = dim(2);
    dAB = dA*dB;
    
    cvx_begin sdp quiet
    variable P(dAB^2, dAB^2) hermitian %AiBiAoBo
    variable k
    
    P_TB = PartialTranspose(P, [2 4], [dA, dB, dA, dB]);
    PAB = PartialTrace(P, [3 4], [dA, dB, dA, dB]);
    
    minimize k
    subject to
        P >= 0;
        -P_TB <= PartialTranspose(JN, [2 4], [dA, dB, dA, dB]) <= P_TB;
        -k*eye(dAB) <= PAB <= k*eye(dAB);
    cvx_end
    l1 = log2(k);
    
    cvx_begin sdp quiet
    variable P(dAB^2, dAB^2) hermitian %AiBiAoBo
    variable k
    
    P_TB = PartialTranspose(P, [2 4], [dA, dB, dA, dB]);
    PAB = PartialTrace(P, [3 4], [dA, dB, dA, dB]);
    
    minimize k
    subject to
        P >= 0;
        -P_TB <= PartialTranspose(JN, [2 4], [dA, dB, dA, dB]) <= P_TB;
        -k * eye(dAB) <= PartialTranspose(PAB, 2, [dA, dB]) <= k * eye(dAB);
    cvx_end

    l2 = log2(k);
    val = max([l1, l2]);
end 