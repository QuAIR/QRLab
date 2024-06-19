function [J1, J2, a] = TCDecompose(JF, dim_list)
    % .. math::
    %
    %     \min\{a:\mathcal{J}_F = a(\mathcal{J}_1 - \mathcal{J}_2)\} 
    %
    % Args:
    %     JF: The input matrix to be decomposed.
    %     dim_list: List specifying the dimensions for the decomposition.
    % 
    % Returns:
    %     [matrix, matrix, numeric]:
    %
    %         J1: The choi matrix of the channel :math:`\mathcal{D}_1`.
    %
    %         J2: The choi matrix of the channel :math:`\mathcal{D}_2`.
    % 
    %         a: The sampling overhead 


    dim = numel(JF(:, 1));

    cvx_begin sdp quiet
    variable J1(dim, dim) complex semidefinite % choi matrix of the positive part of the decoding map
    variable J2(dim, dim) complex semidefinite
    variable a nonnegative

    minimize a
    subject to
        PartialTrace(J1 + J2, 2, dim_list) == a * eye(dim_list(1));
        J1 - J2 == JF;

        % J1(2, 2) == 0;

    cvx_end

    J1 = full(J1);
    J2 = full(J2);
end
