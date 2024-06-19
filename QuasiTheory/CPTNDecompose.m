function [J1, J2, p1, p2] = CPTNDecompose(JF, dim_list)
    % .. math::
    %
    %     \min\{p_1 + p_2:\mathcal{J}_F = p_1\mathcal{J}_1 - p_2\mathcal{J}_2\} 
    %
    % Args:
    %     JF: The input matrix to be decomposed.
    %     dim_list: List specifying the dimensions for the decomposition.
    % 
    % Returns:
    %     [matrix, matrix, numeric, numeric]:
    %
    %         J1: The choi matrix of the channel :math:`\mathcal{D}_1`.
    %
    %         J2: The choi matrix of the channel :math:`\mathcal{D}_2`.
    % 
    %         p1: The sampling overhead of J1
    % 
    %         p2: The sampling overhead of J2

    dim = numel(JF(:, 1));
    
    cvx_begin sdp quiet
    variable J1(dim, dim) complex semidefinite % choi matrix of the positive part of the decoding map
    variable J2(dim, dim) complex semidefinite
    variable p1 nonnegative
    variable p2 nonnegative

    minimize p1 + p2
    subject to

        PartialTrace(J1, 2, dim_list) <= p1 * eye(dim_list(1));
        PartialTrace(J2, 2, dim_list) <= p2 * eye(dim_list(1));
        J1 - J2 == JF;

    cvx_end

    J1 = full(J1 / p1);
    J2 = full(J2 / p2);
end
