function S = QutritSstate()
    % Generate S state
    %
    % Returns:
    %   S (numeric): qudit S state.
    %
    % Note:
    %   Emerson, J. (2014). 
    %   The resource theory of stabilizer computation. 
    %   Bulletin of the American Physical Society, 59.
    
    v1 = [0 1 -1]';
    v1 = v1 / norm(v1);
    S = v1 * v1';
end