function T = QutritTstate()
    % Generate qutrit T state
    %
    % Returns:
    %   T (numeric): qudit T state.
    % Note:
    %   Emerson, J. (2014). 
    %   The resource theory of stabilizer computation. 
    %   Bulletin of the American Physical Society, 59.
    
    zeta = exp(2*pi*i/9);
    vt = [zeta 
          1
          zeta^(-1)];
    vt = vt / norm(vt);
    T = vt * vt';
end