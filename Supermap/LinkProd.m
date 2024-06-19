function Jout = LinkProd(JA, JB, DIM)
    % Link Product of Two Quantum channels, where Aout of JA is linked with
    % the Bin of JB.
    %
    % :Required packages:
    %   `QETLAB  <http://www.qetlab.com/Main_Page>`_
    %
    %
    % .. math::
    %
    %   J_{\mathcal{B} \circ \mathcal{A}} = J_{\mathcal{A}} * J_{\mathcal{B}} = Tr_1[(J_{\mathcal{A}} \otimes I_2)
    %   \cdot (I_0 \otimes J_{\mathcal{B}}^{T_1})],
    %
    % Args:
    %     JA (numeric): The choi matrix of the quantum channel A.
    %     JB (numeric): The choi matrix of the quantum channel B.
    %     DIM (int): The dimensions of channel A and channel B.
    %
    % Returns:
    %     numeric: Resulting Choi matrix of the link product of two Choi matrices JA and JB.
    %
    % Raises:
    %     error: If the dimension of the input state does not match with the pure stabilizer state matrix, an error is raised.
    %
    % :Examples:
    %   .. code-block:: matlab
    %
    %       % Link product of two Choi matrices JA and JB:
    %       Jout = LinkProd(JA, JB, [Ain, Aout, Bin, Bout]);


    
% link Aout and Bin
Ain = DIM(1);
Aout = DIM(2);
Bin = DIM(3);
Bout = DIM(4);

assert(Aout == Bin, 'Output dimension of channel A does not match with the input dimension of channel B');

Link = kron(JA, eye(Bout)) * kron(eye(Ain), PartialTranspose(JB,1,[Bin, Bout]));
Jout = PartialTrace(Link, 2, [Ain, Aout, Bout]);
end