function SE = Stabilizer_Renyi_Entropy(rho,alpha)
    %
    % .. math::
    %
    %     M_\alpha(\psi):=(1-\alpha)^{-1} \log \sum_{P \in \mathcal{P}_n} \Xi_P^\alpha(\psi)-\log d
    % 
    % where :math:`\Xi_P^\alpha(\psi) =d^{-1}\tr^2(P\psi)` and :math: `\mathcal{P}_n` denotes the set of all n-qubit Pauli operators.
    %
    % Args:
    %     rho (numeric): The density matrix of a n-qubit quantum state.
    %     alpha (numeric): alpha parameter.
    %
    % Returns:
    %     numeric: The Stabilizer Renyi Entropy of magic states.
    %
    % Raises:
    %     error: If the input state is not a pure state, an error is raised.
    % 
    % Note:
    % Leone, L., Oliviero, S. F., & Hamma, A. (2022). 
    % Stabilizer rényi entropy. 
    % Physical Review Letters, 128(5), 050402.
    

% Check if the input state is a pure state  
if abs(trace(rho^2) - 1) > 1e-6
    error('Input state is not a pure state. The trace of rho^2 is not equal to 1.');
end

dim = size(rho,1); % rho is the density matrix of a pure state.
n = int8(log2(dim));
Pauli_operators = enumeratePaulis(n);
num_Paulis = size(Pauli_operators,3);

SE = 0;
for i = 1:num_Paulis 
    pro = real(trace(Pauli_operators(:,:,i)*rho))^2/dim;
    SE = SE + pro^alpha;

end
SE  = log2(SE)/(1-alpha) - log2(dim);
end
