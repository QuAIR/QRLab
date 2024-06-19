function Chois_CSPO_cell = RandomCSPO(num_s)

    % Sampling num_s random qubit-qubit channels and sift the qubit-qubit
    % CSPO.
    %
    % Args:
    %   num_s (int): The number of sampling random qubit-qubit channels.
    %
    % Returns:
    %   numeric: cell array containing the sift CSPO.
    %
    % Raises:
    %   error: None.
    %
    % :Examples:
    %   .. code-block:: matlab
    %
    %       Chois_CSPO_cell = RandomCSPO(10000);
    %       % Sampling 10000 random qubit-qubit channels and sift the qubit-qubit
    %       % CSPO.
    % Note:
    %   Wang, X., Wilde, M. M., & Su, Y. (2019). 
    %   Quantifying the magic of quantum channels. 
    %   New Journal of Physics, 21(10), 103002.
    
load('Amat2.mat');
A_mat_2 = A;
num_samples = num_s;
num = 0;
Chois_CSPO_cell = {};
for i = 1:num_samples

    J = RandomSuperoperator(2);
    Js(:,:,i) = J;
    RJs(i) = RoM_Chan(J, A_mat_2);
   
    
    if (abs(RJs(i)-1)<10^(-5)) 
        num = num+1;
      
        Chois_CSPO_cell{end+1} = Js(:,:,i);

    end
  
end

end
