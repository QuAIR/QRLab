function C_R = RobustnessCoherentState(rho)
    % Provide the robustness of a coherent state.
    %
    % .. math::
    %    1 + C_R(\rho) = \min_{\sigma}\{\lambda | \rho \le \lambda \sigma\},
    %
    % where :math:`\sigma` is the incoherent state.
    %
    % Args:
    %     rho (matrix): The density matrix of quantum state.
    %
    % Returns:
    %     C_R (numeric): The robustness of the input state.
    %
    % Note:
    %     Napoli, C., Bromley, T. R., Cianciaruso, M., Piani, M., Johnston, N., & Adesso, G. (2016). 
    %     Robustness of coherence: an operational and observable measure of quantum coherence. 
    %     Physical review letters, 116(15), 150502.

dim = length(rho);
cvx_begin sdp quiet
    variable M(dim,dim) hermitian
    variable s

    minimize s
    subject to
        s >= 0;
        rho <= M;
        trace(M) == 1+s;
        M >= 0;
        M == diag(diag(M));
cvx_end
C_R = s;