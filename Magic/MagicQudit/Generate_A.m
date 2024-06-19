function An = Generate_A(dim, num_copy)

    % 
    % .. math::
    %
    %     T_{\mathbf{u}}=\tau^{-a_1 a_2} Z^{a_1} X^{a_2},
    %     \tau=e^{(d+1) \pi i / d}
    %
    %     A_\mathbf{0}=\frac{1}{d} \sum_{\mathbf{u}} T_{\mathbf{u}}, 
    %     A_{\mathbf{u}}=T_{\mathbf{u}} A_\mathbf{0} T_{\mathbf{u}}^{\dagger}.
    %
    % Generate Heisenberg-Weyl Operators and n-copy Phase-Space Point Operators
    %
    % Args:
    %   dim (numeric): dimension of the operators.
    %   num_copy (numeric): number of copies for phase-space point operators.
    %
    % Returns:
    %   An (numeric): cell array containing the generated operators.
    %
    % Raises:
    %   error: None.
    %
    % :Examples:
    %   .. code-block:: matlab
    %
    %       An = Generate_A(3, 2);
    %       % Generate 2-copy phase-space point operators for dimension 3.
    % Note:
    %   Emerson, J. (2014). 
    %   The resource theory of stabilizer computation. 
    %   Bulletin of the American Physical Society, 59.

    d = dim; % dimension
    dp = 1; % used for the field Z_d for generating Heisenberg-Weyl operators
    X  = GenPauli(1,0,d);
    Z  = GenPauli(0,1,d);
    w  = exp((d+1) * pi * 1i / d);
    
    % 1-copy Heisenberg-Weyl operators T
    T1 = {};
    for a1 = 1:d
        for a2 = 1:d
            T1{end+1} = w^(-(a1-dp) * (a2-dp)) * Z^(a1-dp) * X^(a2-dp);
        end
    end
    
    % n-copy phase-space point operators
    c = 1;
    Tn = T1;
    while c < num_copy
        Ttemp = {};
        for i = 1 : length(Tn)
            for j = 1 : d^2
                Ttemp{end + 1} = Tensor(Tn{i}, T1{j});
            end
        end
        Tn = Ttemp;
    c = c+1;
    end
    
    % n-copy A
    An = {};
    A0 = zeros(d^num_copy);
    for i = 1 : length(Tn)
        A0 = A0 + Tn{i};
    end
    A0 = A0/d^num_copy;
    
    for j=1:length(Tn)
        An{ end + 1} = Tn{j} * A0 * Tn{j}';
    end
end