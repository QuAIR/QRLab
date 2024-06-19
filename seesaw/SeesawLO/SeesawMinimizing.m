function [a, J_out] = SeesawMinimizing(JF, J_in, err, m, opt_sys)
    %   CONVENTION: please read before modify this function!
    %
    %   The data of variables are encoded as follows:
    %
    %   - Suppose qp decomposition terms are arranged by decreasing order of the coefficients
    %       - `m` is then the largest index such that the corresponding coefficient is non-negative.
    %       - `n` is the number of decomposed terms
    %   - `JF` is the target channel.
    %   - `err` is the error tolerance for qp simulation.
    %
    %   - `J_in` is the fixed variable and `J_out` is the opt varible in this iteration.
    %   Depending on `opt_sys`, variable `J_in` and `J_out` could have different meanings:
    %   - If `opt_sys` = 1, then `J_in` = `JA` and `J_out = JB`;
    %   - If `opt_sys` = 2, then `J_in` = `JB` and `J_out = JA`.
    %
    %   - `a` is the coefficient of this optimization, where `a(i)` corresponds to
    %       - the coefficient of i-th term for 1 <= i <= m;
    %       - the negative coefficient of i-th term for m <= i <= n.
    %   - For 1 <= i <= n:
    %       - `JA/JB(:, :, i, 1)` corresponds to the positve part of i-th decomposition term;
    %       - `JA/JB(:, :, i, 2)` corresponds to the negative part of i-th decomposition term.
    %   - `dim_A/B` are the dimensions of system A and B, respectively. 
    %

    %% Process data
    % Permute target channel
    JF = PermuteSystems(JF, [1, 3, 2, 4], [2, 2, 2, 2]);

    % Find number of terms
    n = numel(J_in(1, 1, :));
    assert((1 <= m) && (m <= n));
    
    dim_in = sqrt(numel(J_in(:, 1, 1)));
    assert(floor(dim_in) == dim_in);

    % currently the dimensions for input and output data are the same
    %   we leave interface of dim_out for future usage
    dim_out = dim_in;
    assert(floor(dim_out) == dim_out);

    dim_AB = dim_in * dim_out; % dimension for entire system

    %% Optimize system A
    if opt_sys == 1
        
        dim_A = dim_out;
        dim_B = dim_in;
        JB = J_in;

        % set up cvx program
        cvx_begin sdp quiet
        cvx_precision best

        % create variable state
        variable a(n) nonnegative;
        variable JA(dim_A ^ 2, dim_A ^ 2, n) complex semidefinite;

        % Prepare J_SUM
        J_SUM = 0;

        for i = 1:n
            JA_term = JA(:, :, i);
            JB_term = JB(:, :, i);

            if i <= m
                J_SUM = J_SUM + kron(JA_term, JB_term);
            else
                J_SUM = J_SUM - kron(JA_term, JB_term);
            end

        end

        % objective function
        minimize sum(a)

        subject to

            % eps constraints
            J_SUM - JF <= err * eye(dim_AB ^ 2);
            J_SUM - JF >= - err * eye(dim_AB ^ 2);

            % CPTN constraints
            for i = 1:n
                JA_term_sum = JA(:, :, i);
                PartialTrace(JA_term_sum, 2, [2, 2]) <= a(i) * eye(dim_A);
            end

            sum(a(1:m)) - sum(a((m + 1):n)) >= 1;

        cvx_end
        % cvx_optval

        J_out = JA;
    
    %% Optimize system B
    elseif opt_sys == 2

        dim_A = dim_in;
        dim_B = dim_out;
        JA = J_in;

        % set up cvx program
        cvx_begin sdp quiet
        cvx_precision best

        % create variable state
        variable a(n) nonnegative;
        variable JB(dim_B ^ 2, dim_B ^ 2, n) complex semidefinite;

        % Prepare J_SUM
        J_SUM = 0;

        for i = 1:n
            JA_term = JA(:, :, i);
            JB_term = JB(:, :, i);

            if i <= m
                J_SUM = J_SUM + kron(JA_term, JB_term);
            else
                J_SUM = J_SUM - kron(JA_term, JB_term);
            end

        end

        % objective function
        minimize sum(a)

        subject to

            % eps constraints
            J_SUM - JF <= err * eye(dim_AB ^ 2);
            J_SUM - JF >= - err * eye(dim_AB ^ 2);

            % CPTN constraints
            for i = 1:n
                JB_term_sum = JB(:, :, i);
                PartialTrace(JB_term_sum, 2, [2, 2]) <= a(i) * eye(dim_B);
            end

            sum(a(1:m)) - sum(a((m + 1):n)) >= 1;

        cvx_end
        % cvx_optval

        J_out = JB;

    else
        error('The value of opt_sys %s. is incorrect', opt_sys)
    end

    % process results
    for i = 1:n

        if a(i) > 1e-4
            J_out(:, :, i) = J_out(:, :, i) ./ a(i);
        end
    
    end

    % return results
    a = full(a);
    J_out = full(J_out);
    err = full(err);

end
