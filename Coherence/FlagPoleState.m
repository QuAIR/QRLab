function state = FlagPoleState(dim, p)
% Provide a flag pole state
%
% .. math::
%    \ket{\psi} = \sqrt{p} \ket{0} + \sum_{i=1}^{d-1}
%    \sqrt{\frac{1-p}{d-1}}\ket{i}
%
% Args:
%     dim (numeric): The dimension of the system.
%     p (numeric): The amplitude of the flag.
%
% Returns:
%     state (matrix): The density matrix of the flag state.

ket = [sqrt(p)];
for i = 1:dim-1;
    ket(end + 1) = sqrt((1-p)/(dim-1));
end

state = ket.'*ket;