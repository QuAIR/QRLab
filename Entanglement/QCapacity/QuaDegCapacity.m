function val = QuaDegCapacity(JN, dim)

    % .. math::
    %
    %     Q(\Phi) = \max_{\rho\in\mathcal{D}(\mathcal{H})} I_c(\rho, \Phi)
    %     
    % where :math:`I_c(\rho, \Phi)` is the coherent information of the (degradable) channel :math:`\Phi` for input state :math:`\rho`.
    % 
    % Args:
    %     JN (numeric): The Choi matrix of the channel.
    %     dim (numeric): The dimension of the input and output spaces of the channel.
    %
    % Returns:
    %     numeric: The quantum capacity of the (degradable) channel :math:`\Phi`.

    dA = dim(1); 
    dB = dim(2);
    U = chanconv(JN, 'choi', 'isom'); 
    dE = max(size(U)) / dA;

    Ic = @(rho) quantum_cond_entr(W * applychan(U,rho,'isom', [dA dB]) * W', [ne nf], 2)/log(2); 

    % Quantum capacity = maximum of Ic (for degradable channels) 
    % TODO: generalize to other channels without given degrading map
    cvx_begin sdp quiet
    variable rho(na,na) hermitian 
    maximize (Ic(rho)); 
    subject to
        rho >= 0; trace(rho) == 1; 
    cvx_end
    val = cvx_optval;
end

