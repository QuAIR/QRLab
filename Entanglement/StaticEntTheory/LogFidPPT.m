function val = LogFidPPT(rho, varargin)
    
    % .. math::
    %
    %     E_{\operatorname{PPT}}(\rho_{AB}) = \log\max F(\rho, \sigma) \
    %     s.t. \
    %     \sigma \in \operatorname{PPT}(A:B)
    %
    % Args:
    %   rho (numeric): The density matrix of the bipartite state.
    %   varargin (numeric): The array storing dimensions of subsystems A and B.
    %
    % Returns:
    %   numeric: The PPT-entanglement of :math:`\rho_{AB}`.
    %
    % Raises:
    %   error: If either input/output dimension does not match, an error is raised.
    %
    % Note:
    %   Wang, X., & Wilde, M. M. (2023). 
    %   Exact entanglement cost of quantum states and channels under positive-partial-transpose-preserving operations. 
    %   Physical Review A, 107(1), 012429.


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
    variable X(dAB, dAB) complex
    variable sig(dAB, dAB) hermitian

    t = 0.5*real(trace(X + X'));
    maximize t
    subject to
        [rho, X;
          X', sig] >= 0;
        trace(sig) == 1;
        PartialTranspose(sig, 2, [dA, dB]) >= 0;
        sig >= 0;
    cvx_end

    val = log2(t);
end

