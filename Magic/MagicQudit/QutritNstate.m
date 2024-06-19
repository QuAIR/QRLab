function N = QutritNstate()
    % Generate qutrit N state
    %
    % Returns:
    %   N (numeric): qudit N state.
    % Note:
    %   Emerson, J. (2014). 
    %   The resource theory of stabilizer computation. 
    %   Bulletin of the American Physical Society, 59.
    
    n = [-1;2;-1];
    n_rho = n*n';
    N = n_rho/trace(n_rho);
end