function [c1, c2, J1, J2] = ProbErrorCancel(Noise_channel)
    % Produce an HPTP map (Choi) to inverse the Noise_channel.
    %
    % .. math::
    %    \mathcal{N}^{-1} = c_1 \mathcal{D}_1 + c_2 \mathcal{D}_2.
    %
    % Args:
    %     Noise_channel (matrix): The Choi matrix of the input noise channel.
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
    %     Temme, Kristan, Sergey Bravyi, and Jay M. Gambetta. 
    %     Error mitigation for short-depth quantum circuits.
    %     Physical review letters 119.18 (2017): 180509.



% set the dimension of systems
d_i = sqrt(length(Noise_channel));
d_o = sqrt(length(Noise_channel));

JN = Noise_channel;
JI = d_i*MaxEntangled(d_i,1)*MaxEntangled(d_i,1)'; % choi matrix of the indentity map

% The SDP
cvx_begin sdp quiet
    variable J1(d_i*d_o,d_i*d_o) hermitian
    variable J2(d_i*d_o,d_i*d_o) hermitian
    variable c1
    variable c2

    JD = J1 - J2;  % choi matrix of the decoding map
    JF = PartialTrace(kron(PartialTranspose(JN, 2, [d_i, d_o]), eye(d_i)) * kron(eye(d_i), JD), [2], [d_i,d_o,d_i]);% choi matrix of the final channel

    cost = c1+c2; % cost function (related to the sampling cost)
    minimize cost
    subject to
        J1 >= 0; PartialTrace(J1,2) == c1*eye(d_i);
        J2 >= 0; PartialTrace(J2,2) == c2*eye(d_i);
        JF == JI;
cvx_end

overhead = cost;
JD = JD;
c1 = c1;
c2 = -c2;
J1 = J1/c1;
J2 = J2/c2;