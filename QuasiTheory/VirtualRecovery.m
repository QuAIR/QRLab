function val = VirtualRecovery(rho, dim)

    % .. math::
    %
    %     R^v_{rec}(\rho_{ABC}) = \log\min\{c_1+c_2:(c_1\mathcal{N}_1 - c_2\mathcal{N}_2)(\rho_{AB}) = \rho_{ABC}, c_{1,2}\geq 0, \mathcal{N}_{1,2}\in\operatorname{CPTP}(B,B\otimes C)\}
    %
    % Args:
    %   JN (numeric): The Choi matrix of the bipartite channel.
    %   dim (numeric): The array storing input and output dimensions.
    %
    % Returns:
    %   numeric: The virtual recovery of tripartite state.
    %
    % Raises:
    %   error: If either input/output dimension does not match, an error is raised.
    %
    % Note:
    %   Chen, Y. A., Zhu, C., He, K., Jing, M., & Wang, X. (2023). 
    %   Virtual Quantum Markov Chains. 
    %   arXiv preprint arXiv:2312.02031.

    da = dim(1);
    db = dim(2);
    dc = dim(3);
    
    JI = MaxEntangled(da, 0, 1) * MaxEntangled(da, 0, 1)' * eye(da^2) * da;
    
    cvx_begin sdp quiet
    variable JD1(db*db*dc, db*db*dc) hermitian
    variable JD2(db*db*dc, db*db*dc) hermitian
    variable p1
    variable p2
    JD = JD1 - JD2;
    cost = p1 + p2;
    
    JN = PermuteSystems(kron(JI, JD), [1 3 2 4], [da da db db*dc]);
    sigabc = ApplyMap(PartialTrace(rho, 3, [da, db, dc]), JN);
    
    minimize cost
    subject to
        JD1 >= 0; 
        PartialTrace(JD1, 2, [db db*dc]) == p1 * eye(db);
        JD2 >= 0; 
        PartialTrace(JD2, 2, [db db*dc]) == p2 * eye(db);
        sigabc == rhoabc;
    cvx_end

    val = log2(cost);
end

