clear;
%% Test file for SeesawLO algorithm
% decomposition for CNOT
dA = 2;
dB = 2;

CX = [1 0 0 0;
      0 1 0 0;
      0 0 0 1;
      0 0 1 0];

JN = ChoiMatrix({CX});

%% initialize Seesaw algorithm for LO decomposition
% We have the following default hyperparameters for the algorithm, users
% can input their own settings by adding corresponding arguments in the 
% algorithm function.
% ---- hyperparameters ----
% num_pos = 20;                     number of positive terms
% num_neg = 20;                     number of negative terms
% num_itr = 400;                    total number of opt. iteration
% error_criteria = 1e-7;            error threshold to stop opt.
% cost_convergence_criteria = 5e-7; cost threshold to stop cost minimization
% c = 0.0001;                       initial distortion strength
% -------------------------
% Default optimization
% [a, choi_A, choi_B] = SeesawLO(JN, dA, dB);

% PS. if you want to modify the parameters, just enter your values into 
% the corresponding argument in SeesawLO function.
num_pos = 25;
num_neg = 20;

[a, choi_A, choi_B] = SeesawLO(JN, dA, dB, num_pos, num_neg);

%% output channel reconstruction
J_out = 0;
for j=1:num_pos
    J_out = J_out + a(j)*kron(choi_A(:,:,j), choi_B(:,:,j));
end
for j=1:num_neg
    J_out = J_out - a(num_pos+j)*kron(choi_A(:,:,num_pos+j), choi_B(:,:,num_pos+j));
end
err = TraceNorm(PermuteSystems(J_out, [1,3,2,4], [dA ,dA, dB, dB]) - JN);
fprintf('\n The final simpling cost is %d\n', sum(a));
fprintf('The final trace distance between J_out and JN is %d\n', err);
