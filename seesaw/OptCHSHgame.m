function [rhoAB, A, B] = OptCHSHgame(C, dim)

    % .. math::
    %
    %     \beta_{\operatorname{CHSH}}(S) = \frac{1}{4}\sum_{x,y\in\{0,1\}}
    %     (-1)^{C(a,b,x,y)} \cdot \bra{\psi}A_x \otimes B_y \ket{\psi}
    %
    % Args:
    %   C (numeric): The binary matrix of 4Ds of the given quantum strategy.
    %   dim (numeric): The array storing dimensions of subsystems A and B.
    %
    % Returns:
    %   numeric: The optimal input state to perform the CHSH game.
    %   numeric: The array to store the optimal measurement on system A
    %   numeric: The array to store the optimal measurement on system B. 

    dA = dim(1); dB = dim(2);

    % initialize the input state and the quantum measurement on B for further seesaw optimization
    rhoAB = RandomDensityMatrix(dA*dB,0,1);
    B(:,:,1,1) = RandomDensityMatrix(dB,0,1);
    B(:,:,2,1) = eye(dB) - B(:,:,1,1);
    
    B(:,:,1,2) = RandomDensityMatrix(dB,0,1);
    B(:,:,2,2) = eye(dB) - B(:,:,1,2);

    %% SDP formulation
    for k = 1:10
        %% step 1 optimize A measurement
        cvx_begin sdp quiet
        variable A(dA,dA,2,2) hermitian

        dis = 0;
        for a = 1:2
            for b = 1:2
                for x = 1:2
                    for y = 1:2
                        dis = dis + (-1)^(C(a,b,x,y)) * trace(kron(A(:,:,a,x), B(:,:,b,y))*rhoAB);    
                    end
                end
            end
        end

        dis = real(dis);
        maximize dis
        subject to

            for a = 1:2
                for x = 1:2
                    A(:,:,a,x) >= 0;
                end
            end
    
            A(:,:,1,1) + A(:,:,2,1) == eye(dA);
            A(:,:,1,2) + A(:,:,2,2) == eye(dA);
        
        cvx_end 
        
        %% step 2 optimize B measurement
        cvx_begin sdp quiet
        variable B(dB,dB,2,2) hermitian

        dis = 0;
        for a = 1:2
            for b = 1:2
                for x = 1:2
                    for y = 1:2
                        dis = dis + (-1)^(C(a,b,x,y)) * trace(kron(A(:,:,a,x), B(:,:,b,y))*rhoAB);    
                    end
                end
            end
        end

        dis = real(dis);
        maximize dis
        subject to

            for a = 1:2
                for x = 1:2
                    B(:,:,a,x) >= 0;
                end
            end
    
            B(:,:,1,1) + B(:,:,2,1) == eye(dB);
            B(:,:,1,2) + B(:,:,2,2) == eye(dB);
        
        cvx_end 
                
        %% step 3 optimize shared state
        cvx_begin sdp quiet
        variable rhoAB(dA*dB,dA*dB) hermitian

            dis = 0;
            for a = 1:2
                for b = 1:2
                    for x = 1:2
                        for y = 1:2
                            dis = dis + (-1)^(C(a,b,x,y)) * trace(kron(A(:,:,a,x), B(:,:,b,y))*rhoAB);                    
                        end
                    end
                end
            end

            dis = real(dis);
            maximize dis
            subject to
                rhoAB >= 0;
                trace(rhoAB) == 1;
        cvx_end
    end
    rhoAB = rhoAB; 
    A = A; B = B;
    fprintf('The best CHSH score is %d.\n', dis);
end

