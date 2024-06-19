function val = LogNeg(rho, varargin)

    % .. math::
    %
    %     E_{\operatorname{N}}(\rho_{AB}) = \log\|\rho^{T_{B}}_{AB}\|_1 
    %
    % Args:
    %   rho (numeric): The density matrix of the bipartite state.
    %   varargin (numeric): The array storing dimensions of subsystems A and B.
    %
    % Returns:
    %   numeric: The logarithmic negativity of :math:`\rho_{AB}`.
    %
    % Raises:
    %   error: If either input/output dimension does not match, an error is raised.
    %
    % Note:    
    %   Plenio, M. B. (2005). 
    %   Logarithmic negativity: a full entanglement monotone that is not convex. 
    %   Physical review letters, 95(9), 090503.

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
    variable R(dAB, dAB) hermitian

    rho_TB = PartialTranspose(rho, 2, [dA, dB]);
    
    t = real(trace(rho_TB * R));
    maximize t
    subject to
        -eye(dAB) <= R <= eye(dAB);
    cvx_end

    val = log2(t);
end

