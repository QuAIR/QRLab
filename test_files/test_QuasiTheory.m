%% Test the resource of quasitheory 

% Prepare coefficients
I = [1 0; 0 1];
X = [0 1; 1 0];
Y = [0 -1j; 1j 0];
Z = [1 0; 0 -1];

tol = 10^-8;
noise_level = 0.1;

k0 = [1 0; 0 sqrt(1-noise_level)];
k1 = [0 sqrt(noise_level); 0 0];
J_AD = kron(I, k0) * MaxEntangled(2,1)*MaxEntangled(2,1)'*2 * kron(I, k0') + kron(I, k1) * MaxEntangled(2,1)*MaxEntangled(2,1)'*2 * kron(I, k1');

% Check ProbErrorCancel

abs(ProbErrorCancel(J_AD) - (1+noise_level)/(1-noise_level)) <= tol
abs(ProbErrorCancel(DepolarizingChannel(2, 1-noise_level)) - (1+(1-2/2^2)*noise_level)/(1-noise_level)) <= tol
abs(ProbErrorCancel(DepolarizingChannel(3, 1-noise_level)) - (1+(1-2/3^2)*noise_level)/(1-noise_level)) <= tol
abs(ProbErrorCancel(DepolarizingChannel(4, 1-noise_level)) - (1+(1-2/4^2)*noise_level)/(1-noise_level)) <= tol


% Check ProbErrorCancelO

abs(ProbErrorCancelO(J_AD, X) - 1/(sqrt(1-noise_level))) <= tol
abs(ProbErrorCancelO(DepolarizingChannel(2, 1-noise_level), X) - 1/(1-noise_level)) <= tol
abs(ProbErrorCancelO(DepolarizingChannel(4, 1-noise_level), kron(X,Y)) - 1/(1-noise_level)) <= tol

% Check ProbErrorCancelOS
H = (SwapGenerator(2,2) + SwapGenerator(2,2)')/2;
J_DE2 = PermuteSystems(kron(DepolarizingChannel(2, 1-noise_level), DepolarizingChannel(2, 1-noise_level)), [1,3,2,4]);
J_AD2 = PermuteSystems(kron(J_AD, J_AD), [1,3,2,4]);
abs(ProbErrorCancelOS(J_AD2, H) - 1/(1-noise_level)^2) <= tol
abs(ProbErrorCancelOS(J_DE2, H) - 1/(1-noise_level)^2) <= tol