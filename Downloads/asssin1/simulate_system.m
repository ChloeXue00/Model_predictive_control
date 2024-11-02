function [x_history, u_history] = simulate_system(A, B, Q, R, N, x0)
    % N = 0;
    % while 1
    %     % Increment N
    %     N = N+1;
   N=33;
   % Initialize state and input histories
    x_history = zeros(size(A, 1), 2000);
    u_history = zeros(size(B, 2), 2000);
        x = x0;
 %Start from the maximum N and try to find the smallest N that stabilizes the system
       stable_N = max_N; % This will record the smallest N that resulted in a stable system
       stable = false;

    %  Define Omega and Gamma matrices once, outside the loop
   
   % for i = 1:N
    for N = max_N:-1:1
         Omega = zeros(N*4, 4);
        Gamma = zeros(N*4, N);
       for i = 1:N
         Omega((i-1)*4+1:i*4, :) = A^i;
         for j = 1:i
            Gamma((i-1)*4+1:i*4, j) = A^(i-j)*B;
         end
       end
    Q_hat = blkdiag(kron(eye(N), Q)); % Block diagonal Q
    R_hat = kron(eye(N), R); % Block diagonal R
    K_b= -inv(Gamma'*Q_hat*Gamma + R_hat)*Gamma'*Q_hat*Omega;
    K_0 = K_b(1:1, 1:4);
       % u = -K_0 * x;
   

    stable = false;
    poles = eig(A+(B*K_0));
        for i=1:length(poles)
         if poles(i) >= 1
                stable = false;
         else
             stable_N=N;
             fprintf('The minimum N_3a for stable feedback is N=%.0f\n',N)
             break;
         end
    end    
    for k = 1:2000
        current_horizon = min(N, 2000 - k + 1);
        u = -K_0 * x;
        
        

        
        
        % Apply the first control input
          % Assuming you have computed K0
       %  poles = eig(A+(B*K_0));
       %  for i=1:length(poles)
       %   if poles(i) >= 1
       %          stable = false;
       %   end
       % end    
       %  if stable == true
       %    poles = eig(A+(B*K_0));
       %    fprintf('The minimum N_3a for stable feedback is N=%.0f\n',N)
       %   end
        % Store state and input
        x_history(:, k) = x;
        u_history(:, k) = u;

        % Update state
        x = A * x + B * u;
    end
    if stable
        fprintf('The system became stable at iteration k=%.0f with horizon N=%.0f\n', k, N);
    else
        fprintf('The system did not stabilize within the given time frame.\n');
    end
end
                 
            

 
    
    
    
