function val = LogFidBiNeg(rho, varargin)

    % .. math::
    %
    %     E^{1/2}_{\operatorname{N},2}(\rho_{AB}) = \log\max F(\rho,
    %     \sigma) \ s.t. \
    %     \sigma \in \operatorname{PPT}_2(A:B)
    %
    % Args:
    %   rho (numeric): The density matrix of the bipartite state.
    %   varargin (numeric): The array storing dimensions of subsystems A and B.
    %
    % Returns:
    %   numeric: The logarithmic fidelity of bi-negativity of :math:`\rho_{AB}`.
    %
    % Raises:
    %   error: If either input/output dimension does not match, an error is raised.
    % 
    % Note:
    %   Wang, X., Jing, M., & Zhu, C. (2023). 
    %   Computable and Faithful Lower Bound for Entanglement Cost. 
    %   arXiv preprint arXiv:2311.10649.

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
    variable A(dAB, dAB) hermitian
    variable B(dAB, dAB) hermitian
    variable C(dAB, dAB) hermitian
    variable D(dAB, dAB) hermitian

    t = 0.5*real(trace(X + X'));
    maximize t
    subject to
        [rho, X;
          X', PartialTranspose(C - D, 2, [dA, dB])] >= 0;
        PartialTranspose(A - B, 2, [dA, dB]) == C + D;
        trace(A + B) <= 1;
        A >= 0; B >= 0; 
        C >= 0; D >= 0;
    cvx_end

    val = log2(t);
end

