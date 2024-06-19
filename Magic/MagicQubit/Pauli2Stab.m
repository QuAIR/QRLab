function Stab_Vec = Pauli2Stab(A_mat,n_qubit)
    %
    % :Required packages: 
    %   channel_magic v2.0 https://github.com/jamesrseddon/channel_magic
    % 
    % Convert Pauli representations of pure stabilizers in A_mat to stabilizer matrices:
    %
    % Args:
    %   A_mat: Pauli representations of pure stabilizers from the package channel_magic v2.0.
    %   n_qubit (int): Number of qubits.
    %
    % Returns:
    %   numeric: Stabilizer matrices where each column is a pure stabilizer state.
    %
    % Raises:
    %   error: If the number of qubits does not match with the A_mat file provided, an error is raised.
    %
    % :Examples:
    %   .. code-block:: matlab
    %
    %       A_mat_2 = load('Amat2.mat');
    %       Stab_Vec = Pauli2Stab(A_mat_2, 2);
    %       % Convert Pauli representations of 2-qubit pure stabilizers in A_mat
    %       % to stabilizer matrices.
pauli_array = enumeratePaulis(n_qubit);
d2 = size(A_mat,2)
d1 = size(A_mat,1)

assert(d1 ==4^n_qubit, 'Number of qubits does not match with the A_mat file.');
for i=1:d2
    T_s = zeros(2^n_qubit);
    for j = 1:d1
    T_s = T_s + pauli_array(:,:,j)*A_mat(j,i)
    end
    Stablizer(:,:,i) = T_s/(2^n_qubit);

end
%%%%%%%%%%%%%%%DensityMatrix2Vec
len = size(Stablizer,3)
d = size(Stablizer,1)
for i=1:len
    [D, V] = eig(Stablizer(:,:,i))
Stab_Vec(:,i) = D(:,d);
end
end
