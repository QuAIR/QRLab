function matrix = UnitaryChannel(unitary)
    % UnitaryChannel Provide the Choi matrix of a unitary channel
    % 
    % Args:
    %     unitary (numeric): The unitary matrix.
    %
    % Returns:
    %     numeric: The Choi matrix of the unitary channel.
    
    dim = length(unitary);
    
    max_entangled_state = dim*MaxEntangled(dim) * MaxEntangled(dim)';
    
    matrix = kron(eye(dim), unitary) * max_entangled_state * kron(eye(dim), unitary)';