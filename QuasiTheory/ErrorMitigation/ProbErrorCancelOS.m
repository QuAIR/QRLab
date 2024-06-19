function [overhead, t, JD] = ProbErrorCancelOS(Noise_channel, O, varargin)
    % Provide a CPTP map to inverse the Noise_channel with respect to observable :math:`O`
    % with the method called observable shift.
    %
    % .. math::
    %    \mathcal{N}^\dagger\circ\mathcal{D}^\dagger(O) = \frac{1}{f} (O + tI)
    %
    % Args:
    %     Noise_channel (matrix): The Choi matrix of the noise channel.
    %     O (matrix): Observable.
    %
    % Returns:
    %     [numeric, numeric, matrix]:
    %         overhead: The sampling overhead to retrieve classical information, which equals to :math:`f`.
    %
    %         t: The shifted distance.
    %
    %         JD: The choi matrix of the target channel to implement.
    %
    % Note:
    %     Zhao, B., Jing, M., Zhang, L., Zhao, X., Wang, K., & Wang, X. (2023). 
    %     Retrieving non-linear features from noisy quantum states.
    %     arXiv preprint arXiv:2309.11403.


% set the dimension of systems
d_i = sqrt(length(Noise_channel));
d_o = sqrt(length(Noise_channel));

JN = Noise_channel;

%% The SDP
cvx_begin sdp quiet
    variable JD(d_i*d_o,d_i*d_o) hermitian
    variable p
    variable t

    JF = PartialTrace(kron(PartialTranspose(JN, 2),eye(d_o))*kron(eye(d_i),JD), 2,[d_i,d_o,d_i]);
    cost = p;

    minimize cost
    subject to
        JD >= 0; PartialTrace(JD,2) == p*eye(d_i);
        ApplyMap(O, PermuteSystems(JF.', [2, 1])) == O + t*eye(d_o);

cvx_end
overhead = cost;
JD = JD;
t = t;