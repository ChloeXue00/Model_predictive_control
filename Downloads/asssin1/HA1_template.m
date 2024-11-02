clear;clc;close all;
% Use this template for Home Assignment 1.
% You may add as many variables as you want but DO NOT CHANGE THE NAME OF
% THE VARIABLES WHICH ARE ALREADY IN THE TEMPLATE.
% Also, DO NOT REUSE THE VARIABLES WHICH ARE ALREADY IN THE TEMPLATE FOR
% SOMETHING ELSE.
% The grading will be partially automatized, so it is very important that
% the names of the variables are what the grader script expects to see.
% All the variables whose name must not be changed have a suffix related to
% the number of the question, so just pay extra attention to them.
% Keep also in mind that the variables are expected to have NUMERICAL
% VALUES and not strings of characters.

%% Q1a (discretization of a state-space model)

a = -0.9421;
b = 82.7231;
c = 14.2306;
p = -3.7808;
q = 4.9952;
r = 57.1120;

h = 0.1;

Ac = [0 1 0 0
      b 0 0 a
      0 0 0 1
      q 0 0 p];

Bc = [0; c; 0; r];

C = [1 0 0 0];

D = 0;

% Continuous state-space object
sysc = ss(Ac, Bc, C, D);  

% Discrete Model
%sysd = c2d(sysc, h, 'zoh');  % Convert the system to discrete time using zero-order hold
sysd = expm([Ac, Bc; [0 0 0 0], 0]*h);

A_1a = sysd(1:4,1:4);   %'Matrix A for discrete time model';
B_1a = sysd(1:4,5);     %'Matrix B for discrete time model';
C_1a = C;               %'Matrix C for discrete time model';

eig_1a = eig(A_1a);        %'Eigenvalues of A_1a';

%% Q1b
h=0.1;
t = 0.8*h;  % delay value tao
% Delayed model
sysdl = expm([Ac, Bc; [0 0 0 0], 0]*t);
A_t = sysdl(1:4,1:4);
B_t = sysdl(1:4,5);

sysdlh = expm([Ac, Bc; [0 0 0 0], 0]*(h-t));
A_ht = sysdlh(1:4,1:4);
B_ht= sysdlh(1:4,5);

Aa_1b = [A_ht * A_t, A_ht * B_t; [0,0,0,0],0];  %'Matrix A for delayed discrete time model';
Ba_1b = [B_ht; 1];                      %'Matrix B for delayed discrete time model';
Ca_1b = [C_1a,0] ;                       %'Matrix C for delayed discrete time model';

eig_1b = eig(Aa_1b )  ;                  %'Eigenvalues of Aa_1b';

%% Q2a (Dynamic Programming solution of the LQ problem)
h=0.01;

A = [1.0041,    0.0100,         0,   -0.0000;
     0.8281,    1.0041,         0,   -0.0093;
     0.0002,    0.0000,    1.0000,    0.0098;
     0.0491,    0.0002,         0,    0.9629];

B = [0.0007;
     0.1398;
     0.0028;
     0.5605];

C = [1, 0, 0, 0;
     0, 0, 1, 0]; 


Q = eye(4);
Pf = 10*eye(4);
R = 1;

stable=false;
P=Pf;
N= 33; %guess N
K_b=zeros(1, 4, N);
for k = N-1:-1:0
    % Compute the control gain matrix K at step k
    K_b(:,:,k+1) = -inv(R + B' * P * B) * (B' * P * A);

    % Update the cost-to-go matrix P for the next step
    P = Q + A' * P * A - A' * P * B * inv(R + B' * P * B) * B' * P * A;
end    

% finite horizon problem

while 1
eig_2a = eig(A+B*K_b(:,:,1)); 
poles = abs(eig_2a);
        stable = true;
        for i=1:length(poles)
            if poles(i) >= 1
                stable = false;
            end
        end
 if stable==true
            poles = eig_2a;
            disp(['The minimum N for stable system N=' num2str(N)]);
            disp('The poles are:');
            disp(poles); 
            disp(['The feedback gain K is given by; K = [' num2str(K_b(:).') ']']);

          
            break;
           
 end
end
%eig_2a =

    % 0.5584
    % 0.9236
    % 0.9524
    % 0.9996
N_2a = 33; %'Minimum horizon length to stabilize the system';
K_2a = K_b; %'Stabilizing feedback gain'K = [-46.1566      -5.2278     0.015542     0.403136];

%% Q2b

% infinite horizon problem
[P_inf_2b, K_2b, eig_2b] = dare(A, B, Q, R);
P_inf_2b = 1.0e+04 * ~[4.8730    0.5323   -0.1038  -0.1154
    0.5323    0.0587   -0.0114  -0.0127
   -0.1038   -0.0114    0.0125   0.0027
   -0.1154   -0.0127    0.0027   0.0030 ];  %'Stationary solution of the Riccati equation';

P= zeros(4,4);
P_prev = P;
iteration = 0;
convergenceCriterion = 1e-1;
hasConverged = false;

while ~hasConverged
    % Compute the gain matrix K using current P
    K_2c = -inv(R + B'*P*B)*B'*P*A;
    
    % Update P using the discrete-time Riccati equation
    P_new = Q + A'*P*A - A'*P*B*inv(R + B'*P*B)*B'*P*A;
    
    % Check for convergence
    if norm(P_new - P) <= convergenceCriterion
        hasConverged = true;
    else
        P = P_new;
        iteration = iteration + 1;
    end
    
    % Optional: break if too many iterations to avoid infinite loop
    if iteration > 10000
        warning('Iteration limit reached without convergence.');
        break;
    end
end

disp('Iterations needed for convergence:');
disp(iteration);

% Compare with dare result
[P_inf_DP_2b, ~, ~] = dare(A, B, Q, R);

%iteration = 435;
P_inf_DP_2b =1.0e+04 *[
    4.8730    0.5323   -0.1038  -0.1154
    0.5323    0.0587   -0.0114  -0.0127
   -0.1038   -0.0114    0.0125   0.0027
   -0.1154   -0.0127    0.0027   0.0030];
%'Stationary solution of the Riccati equation obtained through Dynamic Programming)';

%% Q2c
Pf=P_inf_DP_2b;
P=Pf;
N_2c = 1; 
K_2c=zeros(1, 4, N);

while true
    P = Pf;
    stability_achieved = true;
    for k = N_2c-1:-1:0
        K_2c = -inv(R + B'*P*B)*B'*P*A;
        P = Q + K_2c'*R*K_2c + (A+B*K_2c)'*P*(A+B*K_2c);
        if any(abs(eig(A+B*K_2c)) >= 1) % Check stability
            stability_achieved = false;
            break;
        end
    end
    
    if stability_achieved
        break; % Exit loop if stability is achieved
    else
        N_2c = N_2c + 1; % Increase N and try again
    end
end

% Display the new N and K
disp(['The shortest N for stability with Pf = P_inf is: ', num2str(N_2c)]);

N_2c =1; %'Minimum horizon length to stabilize the system';
K_2c =[ -67.7456   -7.5229    0.6939    0.9025]; %'Stabilizing feedback gain';

%% Q3a (batch solution of the LQ problem)

    N = 0;
    while 1
        % Increment N
        N = N+1;
    
        % Construct Q_bar
        Q_bar = [];
        for i=1:N-1
            Q_bar = blkdiag(Q_bar,Q);
        end
        Q_bar = blkdiag(Q_bar,Pf);
    
        
        % Construct R_bar
        R_bar = [];
        for i=1:N
            R_bar = blkdiag(R_bar,R);
        end
        
        
        % Construct omega
        omega = [];
        for i=1:N
            omega = [omega; A^i];
        end
    
 
        % Construct gamma
        gamma = [];
        nolla = zeros(length(B),1);
        lowest_row = [];
        
        for i=1:N
            lowest_row = [lowest_row A^(N-i)*B];
        end
        
        gamma = [lowest_row];
        for i=1:N-1
            above_row = [lowest_row(:,i+1:end)];
            for j=1:i
                above_row = [above_row, nolla];
            end
            gamma = [above_row; gamma];
        end
        
    
        % Determine feedback gain
        K_b = -inv(gamma'*Q_bar*gamma + R_bar)*gamma'*Q_bar*omega;
        
        % K_B = K_b(0)
        K_b = K_b(1,:);
        
        
        eig_3a = eig(A+B*K_b(:,:,1)); 
        poles = abs(eig_3a);
        stable = true;
        for i=1:length(poles)
            if poles(i) >= 1
                stable = false;
            end
        end
        if stable==true
            poles = eig_3a;
            disp(['The minimum N for stable system N=' num2str(N)]);
            disp('The poles are:');
            disp(poles); 
            disp(['The feedback gain K is given by; K = [' num2str(K_b(:).') ']']);

          
            break;
           
        end       
       
   end  
   
% finite horizon problem
eig_3a = poles; 
N_3a =N ; %'Minimum horizon length to stabilize the system';
K0_3a =K_b; %'Stabilizing feedback gain';

% disp(eig_3a);
%disp(N_3a);
%% Q4 (Receding horizon control)
Q = eye(4); % Define Q matrix
Pf = 10*eye(4); % Terminal cost

%4sets together
x0 = [pi/38; 0; 0; 0];
%Initialize the state

[x1, u1] = find_sequence_DP(A,B,Q,  1,Pf,40);
[x2, u2] = find_sequence_DP(A,B,Q,  1,Pf,80);
[x3, u3] = find_sequence_DP(A,B,Q,0.1,Pf,40);
[x4, u4] = find_sequence_DP(A,B,Q,0.1,Pf,80);




figure()
subplot(4,1,1)
hold on; grid on
title('x_1 ball angle')
plot(0:100,x1(1,1:101),'linewidth',2)
plot(0:100,x2(1,1:101),'linewidth',2)
plot(0:100,x3(1,1:101),'linewidth',2)
plot(0:100,x4(1,1:101),'linewidth',2)
legend('R=1,N=40','R=1,N=80','R=0.1,N=40','R=0.1,N=80')
xlabel('timestep k')
ylabel('[rad]')


subplot(4,1,2)
hold on; grid on
title('x_2 ball angular speed')
plot(0:100,x1(2,1:101),'linewidth',2)
plot(0:100,x2(2,1:101),'linewidth',2)
plot(0:100,x3(2,1:101),'linewidth',2)
plot(0:100,x4(2,1:101),'linewidth',2)
legend('R=1,N=40','R=1,N=80','R=0.1,N=40','R=0.1,N=80')
xlabel('timestep k')
ylabel('[rad/s]')


subplot(4,1,3)
hold on; grid on
title('x_3 wheel angle')
plot(0:2000,x1(3,1:2001),'linewidth',2)
plot(0:2000,x2(3,1:2001),'linewidth',2)
plot(0:2000,x3(3,1:2001),'linewidth',2)
plot(0:2000,x4(3,1:2001),'linewidth',2)
legend('R=1,N=40','R=1,N=80','R=0.1,N=40','R=0.1,N=80')
xlabel('timestep k')
ylabel('[rad]')


subplot(4,1,4)
hold on; grid on
title('u input')
plot(0:30,u1(1:31),'linewidth',2)
plot(0:30,u2(1:31),'linewidth',2)
plot(0:30,u3(1:31),'linewidth',2)
plot(0:30,u4(1:31),'linewidth',2)
legend('R=1,N=40','R=1,N=80','R=0.1,N=40','R=0.1,N=80')
xlabel('timestep k')
ylabel('[no unit]')



% plot x1, x2, x3 and u as per the instructions

%% Q5 (constrained receding horizon control)
%constraints const=[x1,x2,x3,x4,u]

A = [1.0041,    0.0100,         0,   -0.0000;
     0.8281,    1.0041,         0,   -0.0093;
     0.0002,    0.0000,    1.0000,    0.0098;
     0.0491,    0.0002,         0,    0.9629];

B = [0.0007;
     0.1398;
     0.0028;
     0.5605];

C = [1, 0, 0, 0;
     0, 0, 1, 0]; 


Q = eye(4);
Pf = 10*eye(4);
R = 1;

stable=false;
P=Pf;
N= 33; %guess N
const = [0 1 0 0 8];
n= size(A, 1); %4,n_states 
 m= size(B, 2); %1 ,number fo_inputs
   x1_max = const(1);
    x2_max = const(2);
    x3_max = const(3);
    x4_max = const(4);
    u_max = const(5);
 
  F      = kron([eye(N); zeros(N)], [0 1 0 0;0 -1 0 0]);
  G      = kron([zeros(N); eye(N)], [1; -1]);
  h      = [x2_max*ones(2*N,1); u_max*ones(2*N,1)];
  %bin      = [x2_max*ones(2*N,1); u_max*ones(2*N,1)];
  t=2000;
  x_vector  = [x0 zeros(n,t)]; 
  u_vector  = zeros(m,t);    
  x = x0;
  
  for iter = 1:t
        
        Qbar = blkdiag(kron(eye(N-1),Q),Pf);
        Rbar = kron(eye(N),R);
        
        H    = blkdiag(Qbar,Rbar);
        f    = []';
        I    = eye(n);

        Aeq1 = kron(eye(N),I)+kron(diag(ones(N-1,1),-1),-A);
         Aeq2 = kron(eye(N),-B);
         Aeq  = [Aeq1 Aeq2];
         beq  = [A*x0;zeros(n*(N-1),1)];
         %Ain = blkdiag(F, G);
        Ain  = [F G];
        bin  = h;
         [Z,~] = CRHC(A,B,N,Q,R,Pf,F,G,h,x,n);

        %[Z,VN,exitflag,~,~] = quadprog(2*H,f,Ain,bin,Aeq,beq);

         options = optimoptions('quadprog', 'Display', 'off');
        [Z,~,exitflag,~,~] = quadprog(2*H, f, Ain, bin, Aeq, beq, [], [], [], options);
        if exitflag ~= 1
             disp(['The exitflag is ' num2str(exitflag) ', which means that something is wrong.'])
                pause
        end
        

        % Get next command and next state
        x     = Z(1:n);
        u     = Z(n*N+1);
        % Save state and input trajectory
        x_vector(:,iter+1) = x;
        u_vector(:,iter)   = u;
        
    
  

   
end




figure()
subplot(4,1,1)
hold on; grid on
title('x1 ball angle')
plot(0:100,x1(1,1:101),'linewidth',2)
plot(0:100,x2(1,1:101),'linewidth',2)
plot(0:100,x3(1,1:101),'linewidth',2)
plot(0:100,x4(1,1:101),'linewidth',2)
legend('R=1,N=40','R=1,N=80','R=0.1,N=40','R=0.1,N=80')
xlabel('timestep k')
ylabel('[rad]')

subplot(4,1,2)
hold on; grid on
title('x2 ball angular speed')
plot(0:100,x1(2,1:101),'linewidth',2)
plot(0:100,x2(2,1:101),'linewidth',2)
plot(0:100,x3(2,1:101),'linewidth',2)
plot(0:100,x4(2,1:101),'linewidth',2)
legend('R=1,N=40','R=1,N=80','R=0.1,N=40','R=0.1,N=80')
xlabel('timestep k')
ylabel('[rad/s]')



subplot(4,1,3)
hold on; grid on
title('x3 wheel angle')
plot(0:2000,x1(3,1:2001),'linewidth',2)
plot(0:2000,x2(3,1:2001),'linewidth',2)
plot(0:2000,x3(3,1:2001),'linewidth',2)
plot(0:2000,x4(3,1:2001),'linewidth',2)
legend('R=1,N=40','R=1,N=80','R=0.1,N=40','R=0.1,N=80')
xlabel('timestep k')
ylabel('[rad]')


subplot(4,1,4)
hold on; grid on
title('u input')
plot(0:30,u1(1:31),'linewidth',2)
plot(0:30,u2(1:31),'linewidth',2)
plot(0:30,u3(1:31),'linewidth',2)
plot(0:30,u4(1:31),'linewidth',2)
legend('R=1,N=40','R=1,N=80','R=0.1,N=40','R=0.1,N=80')
xlabel('timestep k')
ylabel('[no unit]')
% same plot as for Q4 but for the constrained problem