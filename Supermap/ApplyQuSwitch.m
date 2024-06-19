function rho_out = ApplyQuSwitch(qs_k, rho_t, rho_c)  
    % :Required packages: 
    %   `QETLAB  <http://www.qetlab.com/Main_Page>`_
    %
    %
    % Two `n`-qubit quantum channels :math:`\mathcal{N}_1` and :math:`\mathcal{N}_2` have Kraus representations :math:`\{E_i\}_i` and :math:`\{F_j\}_j`
    %
    % .. math::
    %     W_{ij} = \ket{0}\bra{0}_c\otimes E_i^{(2)}F_j^{(1)} + \ket{1}\bra{1}_c\otimes F_j^{(1)}E_i^{(2)}
    %
    % Args:
    %     qs_k (numeric): Kraus operators of the quantum switch channel
    %     rho_t (numeric): Choi matrix of the target quantum state
    %     rho_c (numeric): Choi matrix of the control quantum state
    %
    % Returns:
    %     numeric: Resulting quantum state after applying the quantum switch.
    %
    % :Examples:
    %     .. code-block:: matlab
    %
    %         rho_out = ApplyQuSwitch(qs_k, rho, rho_c);
    %         % Apply Quantum Switch to specific target state and control state.

    

    % Create the input state by taking the Kronecker product
    rho_input = kron(rho_t, rho_c);
    
    % Apply the quantum switch using the provided Kraus operators
    rho_out = ApplyMap(rho_input, qs_k);
end
