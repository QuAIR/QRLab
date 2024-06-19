function val = RainsBound(rho, varargin)

    % .. math::
    %
    %     R(\rho_{AB}) = \min\operatorname{D}(\rho_{AB}\|\sigma_{AB}), \
    %     s.t. \
    %       \sigma_{AB} \in \operatorname{PPT'}(A:B).
    %
    % Args:
    %   rho (numeric): The density matrix of the bipartite state.
    %   varargin (numeric): The array storing dimensions of subsystems A and B.
    %
    % Returns:
    %   numeric: The Rains' bound of :math:`\rho_{AB}`.
    %
    % Raises:
    %   error: If either input/output dimension does not match, an error is raised.
    %
    % Note:
    %   Rains, E. M. (2001). 
    %   A semidefinite program for distillable entanglement. 
    %   IEEE Transactions on Information Theory, 47(7), 2921-2933.

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

    cvx_begin sdp quiet
    variable K(dAB, dAB) hermitian; 
    t1 = (quantum_rel_entr(rho, K)/log(2));
    minimize t1
    subject to
        K >= 0;
        trace(K) <= 1;
        TraceNorm(PartialTranspose(K, 2, [dA, dB])) <= 1;
    cvx_end
    
    val = t1;
end


