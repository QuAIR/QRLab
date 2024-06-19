function [c1, c2, J1, J2] = ProbErrorCancelO(Noise_channel, O)
    % Produce an HPTP map (Choi) to inverse the Noise_channel with respect to observable :math:`O`.
    %
    % .. math::
    %    \mathcal{N}^\dagger\circ\mathcal{D}^\dagger(O) = O,
    %
    % Args:
    %     Noise_channel (matrix): The Choi matrix of the noise channel.
    %     O: Observable.
    %
    % Returns:
    %     [numeric, numeric, matrix, matrix]:
    %         c1: The coefficient of the first component.
    %
    %         c2: The coefficient of the second component.
    %
    %         J1: The choi matrix of the quantum channel :math:`\mathcal{D}_1`.
    %
    %         J2: The choi matrix of the quantum channel :math:`\mathcal{D}_2`.
    %
    % Note:
    %     Zhao, X., Zhao, B., Xia, Z., & Wang, X. (2023). 
    %     Information recoverability of noisy quantum states. 
    %     Quantum, 7, 978.



% set the dimension of systems
d_i = sqrt(length(Noise_channel));
d_o = sqrt(length(Noise_channel));

JN = Noise_channel;
%% The SDP
cvx_begin sdp quiet
    variable J1(d_i*d_o,d_i*d_o) hermitian
    variable J2(d_i*d_o,d_i*d_o) hermitian
    variable c1
    variable c2

    JD = J1 - J2;  % choi matrix of the decoding map
    JF = PartialTrace(kron(PartialTranspose(JN, 2),eye(d_o))*kron(eye(d_i),JD), 2,[d_i,d_o,d_i]);

    cost = c1+c2; % cost function (related to the sampling cost)
    minimize cost
    subject to
        J1 >= 0; PartialTrace(J1,2) == c1*eye(d_i);
        J2 >= 0; PartialTrace(J2,2) == c2*eye(d_i);
        ApplyMap(O, PermuteSystems(JF.', [2, 1])) == O;

cvx_end

%% results
c1 = c1;
c2 = -c2;
J1 = J1/c1;
J2 = J2/c2;