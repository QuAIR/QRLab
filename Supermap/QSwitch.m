function [JoutPlus, JoutMinus] = QSwitch(JN, d) 
    % Quantum Switch Choi Matrices for Control System Measurement when we
    % insert two same quantum channels N, and set the control qubit with :math:`|+\rangle`.
    %
    % :Required packages:
    %   `QETLAB  <http://www.qetlab.com/Main_Page>`_
    %
    %
    % Args:
    %     JN (numeric): Choi matrix of the channel inserted into Quantum
    %      Switch.
    %     d (int): Input dimension of the channel.
    %
    % Returns:
    %     numeric:
    %     JoutPlus  - Choi matrix for the quantum switch when the control system is measured on :math:`|+\rangle`;
    %
    %     JoutMinus - Choi matrix for the quantum switch when the control system
    %     is measured on :math:`|-\rangle`
    %
    % Raises:
    %     error: None.
    %
    % :Examples:
    %   .. code-block:: matlab
    %
    %       [JoutPlus, JoutMinus] = QSwitch(JN, d);
    %       % Compute Choi matrices for a quantum switch with control system in :math:`|+\rangle`
    %       % and :math:`|-\rangle` states.


% 123456, 135 input, 246 output
% 145236
KetI = MaxEntangled(d);
W1 = Tensor(KetI, KetI, KetI);
WPlus = 0.5*(W1 + PermuteSystems(W1, [1 4 5 2 3 6], [d d d d d d]));  
WMinus = 0.5*(W1 - PermuteSystems(W1, [1 4 5 2 3 6], [d d d d d d]));


%Calculate the choi matrix of quantum switch itself when the control system
%is measured on :math:`|+\rangle` with permuting system.
JPlus = d^3*PermuteSystems(WPlus*WPlus', [1 3 5 2 4 6], [d d d d d d]); % Choi matrix for Plus

%Calculate the choi matrix of quantum switch itself when the control system
%is measured on :math:`|-\rangle` with permuting system.
JMinus = d^3*PermuteSystems(WMinus*WMinus', [1 3 5 2 4 6], [d d d d d d]); % Choi matrix for Minus

% Insert two same quantum channels N, calculate the linkproduct of two JNs
% and JPlus
LinkPlus = PartialTranspose(JPlus, [2 3 4 5], [d d d d d d]) * PermuteSystems(Tensor(JN, JN, eye(d^2)), [5 2 4 1 3 6]); % JN23, JN45, I16 -> 135246
JoutPlus = PartialTrace(LinkPlus, [2 3 4 5], [d d d d d d]);

% Insert two same quantum channels N, calculate the linkproduct of two JNs
% and JMinus
LinkMinus = PartialTranspose(JMinus, [2 3 4 5], [d d d d d d]) * PermuteSystems(Tensor(JN, JN, eye(d^2)), [5 2 4 1 3 6]); % JN23, JN45, I16 -> 135246
JoutMinus = PartialTrace(LinkMinus, [2 3 4 5], [d d d d d d]);

end