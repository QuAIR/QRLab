%% SDP for quantum state discrimination
clear;

%% Random density operators rho and sigma
d = 2;
rho = RandomDensityMatrix(d);
sigma = RandomDensityMatrix(d);

% a priori probability p for rho and 1 - p for sigma
priori_prob = rand(1);

% The minimized probabilities for incorrect state discrimination
p_error = 0.5*(1 - TraceNorm(priori_prob*rho - (1 - priori_prob)*sigma))

%% Quantum state discrimination SDP
cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best
variable M(d, d) hermitian

p_error_SDP = trace((eye(d) - M)*priori_prob*rho) + trace(M*(1 - priori_prob)*sigma);
minimize p_error_SDP 

subject to
    M >= 0;
    M <= eye(d);

cvx_end

p_error_SDP


