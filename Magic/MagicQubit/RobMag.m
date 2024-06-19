function [RoM,x] = RobMag(rho,Stab)
    % From https://bartoszregula.me/code/magic
    %
    % .. math::
    %
    %     \mathcal{R}(\rho) = \min_{\textbf{q}}\Big\{ \| \textbf{q} \|_1 \ : 
    %     \rho=\sum_i q_i \ketbra{s_i}{s_i} \ , \ketbra{s_i}{s_i} \in \mathrm{STAB}_n \Big\},
    %
    % Args:
    %     rho (numeric): The density matrix of a n-qubit quantum state.
    %     Stab (numeric): A matrix whose columns are pure stabilizer states.
    %
    % Returns:
    %     numeric: The robustness of magic states.
    %
    % Raises:
    %     error: If the dimension of the input state does not match with the pure stabilizer state matrix, an error is raised.
    %
    % Note:
    % Howard, M., & Campbell, E. (2017). 
    % Application of a resource theory for magic states to fault-tolerant quantum computing. 
    % Physical review letters, 118(9), 090501.
    
dim = size(rho,1);
dim_Stab = size(Stab,1);
assert(dim == dim_Stab, 'Dimension of the input state does not match with the pure stabilizer state');
k = size(Stab,2);

cvx_begin quiet
variable x(k)
X = diag(x);
rho == Stab*X*Stab';
RoM = norm(x,1);
minimize RoM
cvx_end

end