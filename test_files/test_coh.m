%% Test the resource of coherence

% Prepare coefficients
tol = 1e-8;

% Robustness of state 
abs(RobustnessCoherentState(FlagPoleState(2,1/2)) - 1) <= tol
abs(RobustnessCoherentState(FlagPoleState(3,1/3)) - 2) <= tol
abs(RobustnessCoherentState(FlagPoleState(4,1/4)) - 3) <= tol
abs(RobustnessCoherentState(0.5 * FlagPoleState(4,1/4) + 0.5*kron(FlagPoleState(2,0.5), FlagPoleState(2,1))) - 2) <= tol

% Robustness of channel 
theta = pi/4;
U = [cos(theta), -sin(theta);
    sin(theta), cos(theta)];

abs(RobustnessCoherentChannel(UnitaryChannel(U)) - 1) <= tol

abs(RobustnessCoherentChannel(kron(eye(2), FlagPoleState(2, 1/2))) - RobustnessCoherentState(FlagPoleState(2,1/2))) <= tol
abs(RobustnessCoherentChannel(kron(eye(2), FlagPoleState(4, 1/4)), [2,4]) - RobustnessCoherentState(FlagPoleState(4,1/4))) <= tol

% Channel simulation using coherence
JN = RandomSuperoperator(2);
abs(SimulateCoherentChannel(UnitaryChannel(U), FlagPoleState(2, 1/2), 0) - 0) <= tol
abs(SimulateCoherentChannel(JN, FlagPoleState(2, 1/(1+RobustnessCoherentChannel(JN))^2), 0) - 0) <= tol

% Probabilistic channel simulation using coherence
theta = pi/8;
U = [cos(theta), -sin(theta);
    sin(theta), cos(theta)];

J_U = PermuteSystems(kron(UnitaryChannel(U), UnitaryChannel(U)), [1,3,2,4]);
abs(SimulateCoherentChannelProb(J_U, FlagPoleState(2,1/2), 0, 0) - 1/((1+sin(2*theta))^2-1)) <= tol
abs(SimulateCoherentChannelProb(RandomSuperoperator(2), FlagPoleState(2,1/2), 0, 0) - 1) <= tol
