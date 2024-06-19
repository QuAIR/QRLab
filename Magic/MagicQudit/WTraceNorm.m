function [output] = WTraceNorm(PhaseOp, rho)
    % Compute the Wigner Trace Norm of a linear Operator
    %
    % Args:
    %   PhaseOp (numeric): Cell array of Phase space operators.
    %   rho (numeric): Density matrix.
    %
    % Returns:
    %   output (numeric): Wigner trace norm of the given operator.
    %
    % Raises:
    %   error: If either input/output dimensions does not match, an error
    %    is raised.
    %
    % :Examples:
    %   .. code-block:: matlab
    %
    %       [output] = WTraceNorm(PhaseOp, rho);
    %       % Calculate the Wigner Trace Norm of rho.
    % 
    % Note:
    %   Wang, X., Wilde, M. M., & Su, Y. (2020). 
    %   Efficiently computable bounds for magic state distillation. 
    %   Physical review letters, 124(9), 090505.

    dim = size(rho);    %dimension
    output=0;
    for j = 1:length(PhaseOp)
        w = abs(real( (trace(PhaseOp{j}*rho) / dim(1) )));
        output = w + output;
    end
end