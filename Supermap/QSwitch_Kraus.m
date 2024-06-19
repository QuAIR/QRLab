function QS_K = QSwitch_Kraus(Kraus_o1,Kraus_o2)
    % Two :math:`n`-qubit quantum channels :math:`\mathcal{N}_1` and :math:`\mathcal{N}_2` have Kraus representations :math:`\{E_i\}_i` and :math:`\{F_j\}_j`
    %
    % .. math::
    %    W_{ij} = \ket{0}\bra{0}_c\otimes E_i^{(2)}F_j^{(1)} + \ket{1}\bra{1}_c\otimes F_j^{(1)}E_i^{(2)}
    %
    % Args:
    %     Kraus_o1 (numeric): Cell array of Kraus operators for the first
    %      channel.
    %     Kraus_o2 (numeric): Cell array of Kraus operators for the second
    %      channel.
    %
    % Returns:
    %     numeric: Cell array of quantum switch Kraus operators.
    %
    % Raises:
    %     error: None.
    %
    % :Examples:
    %   .. code-block:: matlab
    %
    %       QS_K = QSwitch_Kraus(Kraus_o1, Kraus_o2);
    %       % Generate quantum switch Kraus operators from two sets of Kraus
    %       % operators.


    rho_0 = [1. ,0.;0., 0.];
    rho_1 = [0., 0.;0., 1.];
    QS_K = {};
    k=1;
    for i = 1:numel(Kraus_o2)
        for j = 1:numel(Kraus_o1)
            superoperator = kron(Kraus_o2{i} * Kraus_o1{j}, rho_0) + kron(Kraus_o1{j} * Kraus_o2{i}, rho_1);
            QS_K{k} = superoperator; % Append to cell array
            k = k+1;
        end
    end
end