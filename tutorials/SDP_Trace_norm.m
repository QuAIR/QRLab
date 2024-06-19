%% SDP for trace norm
clear;

%% Random Hermitian operator H
d = 2;
a = rand(2,2) + rand(2,2)*1i;
H = a*a';

% Trace norm
trace_norm = TraceNorm(H)

%% Trace norm primal SDP
cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best
variable lambda1(d, d) hermitian
variable lambda2(d, d) hermitian

trace_norm_primal_SDP = trace(H*(lambda1 - lambda2));
maximize trace_norm_primal_SDP 

subject to
    lambda1 >= 0;
    lambda2 >= 0;
    lambda1 <= eye(d);
    lambda2 <= eye(d);

cvx_end

trace_norm_primal_SDP

%% Trace norm dual SDP
cvx_begin sdp quiet
cvx_solver sedumi
cvx_precision best
variable Y1(d, d) hermitian
variable Y2(d, d) hermitian

trace_norm_dual_SDP = trace(Y1 + Y2);
minimize trace_norm_dual_SDP 

subject to
    Y1 >= 0;
    Y2 >= 0;
    Y1 >= H;
    Y2 >= -H;

cvx_end

trace_norm_dual_SDP
