function mana = Chan_Mana(JN, DIM)

    % .. math::
    %
    %     \mathcal{M}(\mathcal{N}_{A \rightarrow B}) = 
    %     \log \max_\mathbf{u} \| \mathcal{N}_{A \rightarrow B}(A_A^\mathbf{u})\|_{W,1} 
    %     = \log \max_u \sum_v |W_{\mathcal{N}}(\mathbf{v}|\mathbf{u})|
    %
    % Args:
    %   JN (numeric): The choi matrix of the quantum channel.
    %   DIM (int): The dimension vector [da db], where da and db are the input and
    %                   output dimensions of the channel.
    %
    % Returns:
    %   numeric: Channel's mana.
    %
    % Raises:
    %   error: If the input and output dimension does not match with the channel, an error is raised.
    %
    % :Examples:
    %   .. code-block:: matlab
    %
    %       mana = Chan_Mana(JN, [3 3]);
    %       % Compute the mana of a qutrit-qutrit channel.
    % 
    % Note:
    %   Wang, X., Wilde, M. M., & Su, Y. (2019). 
    %   Quantifying the magic of quantum channels. 
    %   New Journal of Physics, 21(10), 103002.

    da = DIM(1);
    db = DIM(2);
    % Generate the set of all possible phase space operators of dimension da
    A_in = Generate_A(da, 1);
    A_out = Generate_A(db, 1);
    mana = 0;
    assert(da*db == size(JN,1), 'Dimension does not match with the channel');
    for u = 1:length(A_in)
        m = WTraceNorm(A_out, ApplyMap(A_in{u}, JN));
        if m >= mana
            mana = m;
        end
    end
    mana = log2(mana);
end
