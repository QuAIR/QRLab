function val = EntAssistedCapacity(JN, dim)

    % .. math::
    %
    %     C_{ea}(\Phi) = \max_{\rho\in\mathcal{D}(\mathcal{H})} I(\rho, \Phi)
    %     
    % where :math:`I(\rho, \Phi)` is the mutual information of the channel :math:`\Phi`.
    % 
    % Args:
    %     JN (numeric): The Choi matrix of the channel.
    %     dim (numeric): The dimension of the input and output spaces of the channel.
    %
    % Returns:
    %     numeric: The entanglement-assisted classical capacity of :math:`\Phi`.
    %
    % Note:
    %     Bennett, C. H., Shor, P. W., Smolin, J. A., & Thapliyal, A. V. (2002). 
    %     Entanglement-assisted capacity of a quantum channel and the reverse Shannon theorem. 
    %     IEEE transactions on Information Theory, 48(10), 2637-2655.

    % Dimensions of input, output, and environment spaces of channel 
    dA = dim(1); 
    dB = dim(2);
    U = chanconv(JN, 'choi', 'isom'); 
    dE = max(size(U)) / dA;

    cvx_begin sdp quiet
    variable rho(dA, dB) hermitian; 
    Urho = U * rho * U';
    maximize quantum_cond_entr(Urho, [dB dE]) + quantum_entr(PartialTrace(Urho, 2, [dB dE]))/log(2); 
    subject to
        rho >= 0; trace(rho) == 1; 
    cvx_end
    val = cvx_optval;
end

