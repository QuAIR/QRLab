function val = MaxRains(rho, varargin)
    
    % .. math::
    %
    %     E_{\operatorname{W}}(\rho_{AB}) =
    %     \log\max\operatorname{Tr}[\rho^{T_B}_{AB} R_{AB}], \ s.t.\
    %       \|R_{AB}\|_{\infty} \leq 1, R^{T_{B}}_{AB} \geq 0.
    %
    % Args:
    %   rho (numeric): The density matrix of the bipartite state.
    %   varargin (numeric): The array storing dimensions of subsystems A and B.
    %
    % Returns:
    %   numeric: The improved logarithmic negativity (max-Rains entropy) of :math:`\rho_{AB}`.
    %
    % Raises:
    %   error: If either input/output dimension does not match, an error is raised.
    %
    % Note:
    %   Wang, X., & Duan, R. (2016). 
    %   Improved semidefinite programming upper bound on distillable entanglement. 
    %   Physical Review A, 94(5), 050301.

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
        PartialTranspose(R, 2, [dA, dB]) >= 0;
    cvx_end

    val = log2(t);
end

