function H = QutritHstate()
   % Generate H state
    %
    % Returns:
    %   H (numeric): qudit H state.
    %
    % Note:
    %   Emerson, J. (2014). 
    %   The resource theory of stabilizer computation. 
    %   Bulletin of the American Physical Society, 59.
    
    wh = exp(2*pi*i/3);
    H_gate = 1/sqrt(3) * [1 1 1; 1 wh wh^2; 1 wh^2 wh];
    [S D] = eig(H_gate);
    v = S(:,1);
    H = v*v';
end