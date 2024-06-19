function [cvx_optval, X] = TempLogNeg(rho, varargin)
    
    % .. math::
    %
    %     E_N^{\tau}(\rho_{AB}).
    %
    % Args:
    %   rho (numeric): The density matrix of the bipartite state.
    %   varargin (numeric): The array storing dimensions of subsystems A and B.
    %
    % Returns:
    %   numeric: The Temperal logarithmic negativity of :math:`\rho_{AB}`.
    %
    % Raises:
    %   error: If either input/output dimension does not match, an error is raised.
    %
    % Note:
    %   Lami, L., & Regula, B. (2023). 
    %   No second law of entanglement manipulation after all. 
    %   Nature Physics, 19(2), 184-189.

    if nargin<2
        omega = rho;
        d = sqrt(max(size(rho)))*[1,1];
    elseif nargin==2
        if min(size(varargin{1})) == 1
            omega = rho;
            d = varargin{1} * ones(1,3-max(size(varargin{1})));
        else
            omega = varargin{1};
            d = sqrt(max(size(rho)))*[1,1];
        end
    elseif nargin==3
        omega = varargin{1};
        d = varargin{2} * ones(1,3-max(size(varargin{2})));
    end
    
    dd = prod(d);
    
    rho = (rho + rho') / 2; % to avoid numerical issues
    omega = (omega + omega') / 2;
    
    cvx_begin sdp quiet
    variable X(dd,dd) hermitian
    maximize trace(X*rho)
    subject to
        -eye(dd) <= PartialTranspose(X, 2, d) <= eye(dd);
        -trace(X * omega) * eye(dd) <= X <= trace(X * omega) * eye(dd);
    cvx_end

    cvx_optval = log2(cvx_optval);
end

