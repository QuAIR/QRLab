function R_star = RoM_Chan(J,A_mat)
    % :Dependencies:
    %   Trans_K
    %   findChoiCPTProbustness  from channel_magic v2.0 https://github.com/jamesrseddon/channel_magic
    % 
    % .. math::
    %
    %     \mathcal{R}_*(\mathcal{N}) := \min_{\mathcal{N}_{\pm} \in \text{CSPO}} \big\{ 2p+1: \mathcal{N}=(1+p)\mathcal{N}_+ - p \mathcal{N}_-, p\geq 0\big\}
    %
    % Compute the Channel robustness of a Choi matrix
    %
    % Args:
    %   J (numeric): Choi matrix of a quantum channel.
    %   A_mat: Pauli representations of pure stabilizers form the package channel_magic v2.0.
    %
    % Returns:
    %   numeric: Channel robustness of the channel.
    %
    % Raises:
    %   error: If dimension of channel does not match with the A_mat file provided, an error is raised.
    %
    % Note:
    % Seddon, J. R., & Campbell, E. T. (2019). 
    % Quantifying magic for multi-qubit operations. 
    % Proceedings of the Royal Society A, 475(2227), 20190251.

    
    

dim_J = size(J,1);
dim_A = size(A_mat,1);
assert(dim_J == dim_A, 'Dimension of the channel does not match with the A_mat file provided');
JK = Trans_K(J);
[R_star,~,distribution_JK,~,~] = ...
    findChoiCPTProbustness(JK,A_mat,'high');
end