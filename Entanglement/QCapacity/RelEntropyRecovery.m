function [val, JR] = RelEntropyRecovery(rhoABC, dim)

    % .. math::
    %
    %     R_{rec}(\rho_{ABC}) = \min_{\mathcal{R}:B\rightarrow BC} D(\rho_{ABC}\| (\mathcal{I}_A \otimes \mathcal{R})(\rho_{AB}))
    % 
    % Args:
    %     rhoABC (numeric): The density matrix of the tripartite target state.
    %     dim (numeric): The dimension of the composite systems of the state.
    %
    % Returns:
    %   [numeric, numeric]:
    %     val: The relative entropy of recovery of :math:`\rho_{ABC}`.
    % 
    %     JR: The Choi matrix of the recovery channel.

    assert(numel(dim) >= 3, 'The number of subsystems must be 3!');

    dA = dim(1);
    dB = dim(2);
    dC = dim(3);
    rho_dim = size(rhoABC);
    assert(dA*dB*dC == rho_dim(1), 'State dimension does not match!');

    JI = MaxEntangled(dA, 0, 1) * MaxEntangled(dA, 0, 1)' * eye(dA^2) * dA;

    rhoAB = PartialTrace(rhoABC, 3, [dA dB dC]); 

    cvx_begin sdp quiet
    variable JR_B_BC(dB^2*dC,dB^2*dC) hermitian 
    variable chanout_ABC(dA*dB*dC, dA*dB*dC) hermitian 

    minimize quantum_rel_entr(rhoABC, chanout_ABC) / log(2)
    subject to

        JR_B_BC >= 0; 
        PartialTrace(JR_B_BC, 2, [dB dB*dC]) == eye(dB); 

        JR_AB_ABC = PermuteSystems(kron(JI, JR_B_BC), [1 3 2 4], [dA dA dB dB*dC]);
        chanout_ABC = applychan(JR_AB_ABC, rhoAB, 'choi2', [dA*dB dA*dB*dC]);
    cvx_end

    val = cvx_optval;
    JR = JR_B_BC;
end

