function H = QubitHstate()
   % Generate H state
    %
    % Returns:
    %   H (numeric): qubit H state.
    %
    % Note:
    %   Emerson, J. (2014). 
    %   The resource theory of stabilizer computation. 
    %   Bulletin of the American Physical Society, 59.
    
    h = [cos(pi/8); sin(pi/8)];
    H = h*h';
end