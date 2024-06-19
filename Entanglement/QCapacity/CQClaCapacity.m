function val = CQCapacity(map_alph)

    % .. math::
    %
    %     C(\Phi) = \max S\left(\sum_{x\in\mathcal{X}} p_x \Phi(x)\right) - \sum_{x\in\mathcal{X}} p_x S(\Phi(x))
    %     
    % where :math:`{p_x}_x` forms probability distribution, :math:`S` is the von Neumann entropy, 
    % and :math:`\Phi:\mathcal{X}\rightarrow \mathcal{D}(\mathcal{H})` is some classical-quantum channel.
    % 
    % Args:
    %     map_alph (numeric): mapping data of the CQ-channel, should be stored in a 3D-tensor with the 
    %     first index representing the labelling :math:`x\in \mathcal{X}` and last two indices give 
    %     mapped density operators.
    %
    % Returns:
    %     numeric: The capacity value of :math:`\Phi`.
    %
    % Note:
    %     Hayashi, M., & Nagaoka, H. (2003). 
    %     General formulas for capacity of classical-quantum channels. 
    %     IEEE Transactions on Information Theory, 49(7), 1753-1768.

    map_size = size(map_alph);
    num_density = map_size(1);

    cvx_begin sdp quiet
        variables p(num_density) 

        compose_state = zeros(map_size(2:3));
        indiv_entropy = 0;
        for x = 1 : num_density
            compose_state = compose_state + p(x) * map_alph(x,:,:);
            indiv_entropy = indiv_entropy + p(x) * quantum_entr(map_alph(x,:,:))/log(2);
        end
        compose_entropy = quantum_entr(compose_state)/log(2);

        maximize compose_entropy - indiv_entropy
        subject to
            sum(p) == 1;
            p >= 0; 
    cvx_end
    val = compose_entropy - indiv_entropy;

end

