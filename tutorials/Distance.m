%% Trace Distance and Fidelity
clear;

% Intialize a quantum density operator
v0 = [1,0];
d_in = v0'*v0;

% Define a depolarizing channel with probability 0.5: 
% Depol(X) = (1-p)*Tr(X)*I/d^2 + p*X
Depol = DepolarizingChannel(2, 0.5);

% Apply the channel to the state
d_out = ApplyMap(d_in, Depol);

% Trace distance between d_out and d_in
dis = 0.5*TraceNorm(d_out - d_in);
dis

% Fidelity between d_out and d_in
fid = Fidelity(d_in, d_out);
fid

%% Diamond Norm

% The Choi matrix for identity channel
Iden = ChoiMatrix({eye(2), eye(2)});

% Diamond norm between p = 0.5 depolarizing channel and identity channel
prob = 0.5;
Depol = DepolarizingChannel(2, prob);
diamond_norm = DiamondNorm(Depol - Iden);
diamond_norm


% The diamond norm between depolarizing channel and identity channel increases as
% prob increases
prob = 0:0.05:1;
diamond_norm = [];
for i = 1:1:length(prob)
Depol = DepolarizingChannel(2, 1 - prob(i));
diamond_norm(i) = DiamondNorm(Depol - Iden);
end

plot(prob, diamond_norm)
xlabel('prob')
ylabel('Diamond norm')
title('Diamond norm between depolarizing channel and identity channel')