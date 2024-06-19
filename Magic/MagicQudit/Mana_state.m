function mana = Mana_state(rho)
    % .. math::
    %
    %     \mathcal{M}(\rho) = \log \sum_\mathbf{u} |W_{\rho}(\mathbf{u})|,
    %
    % where :math:`W_{(\rho)}` is the Wigner representation of :math:`\rho`. 
    %
    % Args:
    %   rho (numeric): The density matrix of the quantum state.
    %
    % Returns:
    %   numeric: Magic mana of the given quantum state.
    %
    % Raises:
    %   error: If either input/output dimension does not match, an error is raised.
    % 
    % Note:
    %   Emerson, J. (2014). 
    %   The resource theory of stabilizer computation. 
    %   Bulletin of the American Physical Society, 59.
    
    % Calculate the mana of a given state
    DIM = size(rho);
    d = DIM(1);
    A = Generate_A(d, 1);
    
    % calculate magic mana of the given state
    SUM = 0;
    for i=1:length(A)
        w = trace(A{i}*rho)/d;
        SUM = SUM + abs(real(w));
    end

    mana = log2(SUM);
end

