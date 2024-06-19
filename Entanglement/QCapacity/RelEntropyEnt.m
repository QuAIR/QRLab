function val = RelEntropyEnt(rho, dim)

    % .. math::
    %
    %     E_{\operatorname{PPT}}(\rho) = \min_{\tau\in\operatorname{PPT}} D(\rho\| \tau)
    % 
    % Args:
    %     rho (numeric): The density matrix of the target state.
    %     dim (numeric): The dimension of the composite systems of the state.
    %
    % Returns:
    %     numeric: The PPT relative entropy of :math:`\rho`.

    dA = dim(1); 
    dB = dim(2);
    rho_dim = size(rho);
    assert(dA*dB == rho_dim(1), 'State dimension does not match!');

    cvx_begin sdp quiet
    variable tau(dA*dB, dA*dB) hermitian 
    minimize quantum_rel_entr(rho, tau)/log(2)
    subject to 
        tau >= 0; trace(tau) == 1; 
        PartialTranspose(tau, 2, [dA dB]) >= 0; % Positive partial transpose constraint 
    cvx_end
    val = cvx_optval;
end

