function val = MaxRainsInfo(JN, dim)
    
    % .. math::
    %
    %     R^{2\rightarrow 2}_{\max}(\mathcal{N}) = \log \inf\{\|V_{AB} + Y_{AB}\|_{\infty}: 
    %      (V_{ABA'B'} - Y_{ABA'B'})^{T_{BB'}} \geq J^{\mathcal{N}}_{ABA'B'}\},
    %
    % Args:
    %   JN (numeric): The Choi matrix of the bipartite channel.
    %   dim (numeric): The array storing input and output dimensions.
    %
    % Returns:
    %   numeric: The bidirectional max-Rains information of bipartite channel.
    %
    % Raises:
    %   error: If either input/output dimension does not match, an error is raised.
    %
    % Note:
    %   Wang, X., Fang, K., & Duan, R. (2018). 
    %   Semidefinite programming converse bounds for quantum communication. 
    %   IEEE Transactions on Information Theory, 65(4), 2583-2592.

    dA = dim(1);
    dB = dim(2);
    dAB = dA*dB;
    
    cvx_begin sdp quiet
    variable V(dAB*dAB,dAB*dAB) hermitian
    variable Y(dAB*dAB,dAB*dAB) hermitian
    variable k

    VAB = PartialTrace(V, [3,4], [dA, dB, dA, dB]);
    YAB = PartialTrace(Y, [3,4], [dA, dB, dA, dB]);

    minimize k
    subject to
        V >= 0; Y >= 0; 
        PartialTranspose(V - Y, [2,4], [dA, dB, dA, dB]) >= JN;
        -k * eye(dAB) <= VAB + YAB <= k * eye(dAB);
    cvx_end 
    val = log2(k);
end