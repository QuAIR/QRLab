function out = Eeta(rho, varargin)

    % .. math::
    %
    %     E_{\eta}(\rho_{AB}) = \max - \log\|Y_{AB}^{T_{B}}\|_{\infty} \
    %     s.t. \
    %     -Y_{AB} \leq P^{T_{B}} \leq Y_{AB}
    % 
    % where :math:`P` is the projection onto :math:`\operatorname{supp}(\rho_{AB})`. 
    %
    % Args:
    %     rho (numeric): The density matrix of the bipartite state.
    %     varargin (numeric): The array storing dimensions of subsystems A and B.
    %
    % Returns:
    %     numeric: :math:`\eta`-entanglement of :math:`\rho_{AB}`.
    %
    % Raises:
    %     error: If either input/output dimension does not match, an error is raised.
    %
    % Note:
    %     Wang, X., & Duan, R. (2017). 
    %     Irreversibility of asymptotic entanglement manipulation under quantum operations completely preserving positivity of partial transpose. 
    %     Physical Review Letters, 119(18), 180506.

    if nargin<2
        omega = rho;
        d = sqrt(max(size(rho)))*[1,1];
    elseif nargin==2
        if min(size(varargin{1})) == 1
            omega = rho;
            d = varargin{1} * ones(1,3-max(size(varargin{1})));
        else
            omega = varargin{1};
            d = sqrt(max(size(rho)))*[1,1];
        end
    end

    dA = d(1);
    dB = d(2);
    dAB = dA*dB;

    [V, D] = eig(rho);
    P = 0;
    for i = 1:dAB
        if D(i, i) > 0.00001
            P = P + V(:, i) * V(:, i)';
        end
    end

    cvx_begin sdp quiet
    variable Y(dAB, dAB) hermitian
    
    YB = PartialTranspose(YB, 2, [dA dB]);
    f = norm(YB, inf);

    minimize f
    subject to
        PartialTranspose(P, 2, [dA dB]) + Y >= 0;
        PartialTranspose(P, 2, [dA dB]) - Y <= 0;
    cvx_end
    
    out = -log2(f);
end