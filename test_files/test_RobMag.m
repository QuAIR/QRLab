% clear;
load('Amat1.mat'); % Load A matrix for all 1-qubit stabiliser states.
A_mat_1 = Amat1;


d = 2;

% JT1 = Chois_CSPO(:,:,5); 
% JT2 = Chois_CSPO(:,:,5);

rho = [1,0;0,0];


Stab = Pauli2Stab(A_mat_1,1);
[R,x] = RobMag(rho,Stab);
