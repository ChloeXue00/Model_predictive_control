clc
close all
clear

%% Q1a
% find terminal weight to guarantee asymptotic stability and plot X_N


% Autonomous system
% x(k+1) = A*x(k)+ B*u(k)
A=[1.2 1;0 1]; % asymptotically stable
B=[0;1];

% create system in MPT3 interface
sys = LTISystem('A',A,'B',B);

% constraints on states
sys.x.min = [-15;-15];
sys.x.max = [ 15; 15];
sys.u.min = -1;
sys.u.max = 1;
X = Polyhedron('lb',sys.x.min,'ub',sys.x.max);
U = Polyhedron('lb',sys.u.min,'ub',sys.u.max);
Q=eye(2);
R=100;
N=4;

% Set model penalties
sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);

% Xf=0 as terminal set
Xf = Polyhedron([0 0]);
sys.x.with('terminalSet');
sys.x.terminalSet = Xf;

[Pf, L, G] = dare(A, B, Q, R);
sys.x.with('terminalPenalty');
sys.x.terminalPenalty = QuadFunction(Pf);
% Compute LQR gain K
[K, S, e] = dlqr(A, B, Q, R);

% Check for stability
eig_closed_loop = eig(A - B*K);

% Display eigenvalues of closed-loop system
disp('Eigenvalues of the closed-loop system:');
disp(eig_closed_loop);

%Check if all eigenvalues are inside the unit circle for stability
is_stable = all(abs(eig_closed_loop) < 1);

if is_stable
    disp('The system is asymptotically stable.');
else
    disp('The system is NOT asymptotically stable.');
end

figure(1);
% XN = Xf;
for i = 1:N
    X = sys.reachableSet('X', Xf, 'U', U, 'N',i,'direction', 'backward');
    XN = X.intersect(X).minHRep();
    
    XN.plot('color','blue','alpha',(0.1*i));
    title('XN');
    hold on;
end
xlabel('$x_1$')
ylabel('$x_2$')
legend('X_N, N=1','X_N, N=2','X_N, N=3','X_N, N=4','(0,0)')
title('X_N')


disp('Pf=');
disp(Pf);



Pf_1a =Pf;% 'terminal weigth to guarantee asymptotic stability';
X_N_1a =XN;% 'set of feasible initial states for horizon length N (to be plotted)';

%% Q1b
% design 3 different MPC controllers, simulate them and plot the results
% (both for the simulated states and for the actual ones)


x0=[7;-4];

N1=10;
N2=15;
N3=20;
N  = [N1,N2,N3]; % RHC horizon length 


%RHC reachable states:

for i = 1:N1 
    % compute backward reachable set
    XN1=sys.reachableSet('X',XN,'U',U,'direction','backward');
    % intersect with the state constraints
    XN1 = XN1.intersect(X).minHRep(); 
end


for i = 1:N2 
  
    XN2 = sys.reachableSet('X',XN,'U',U,'direction','backward');
    XN2 = XN2.intersect(X).minHRep();
    
end


for i = 1:N3 
    XN3 = sys.reachableSet('X',XN,'U',U,'direction','backward');
    XN3 = XN3.intersect(X).minHRep();
  
end

% MPC controller simulation,state trajectories
%mpc1
mpc1=MPCController(sys,N1);
loop1 = ClosedLoop(mpc1,sys);
data1=loop1.simulate(x0,N1);
 [~,feasible,openloop1]=mpc1.evaluate(x0);
 [~,feasible,closedloop1]=mpc1.evaluate(x0);
    if feasible==0
        disp('not feasible')
        
    end
figure(6)
subplot(2,1,1);
stairs(0:N1-1,data1.U,'MarkerSize',12,'LineWidth',1);
hold on;grid on;
plot(data1.X(1,:),'MarkerSize',12,'LineWidth',1);
hold on;
plot(data1.X(2,:),'MarkerSize',6,'LineWidth',1);
xlabel('iteration');
% legend('x_1','x_2','u');
% legend('MPC x_1','MPC x_2','MPC u');
title('MPC1 state trajectories VS prediction');
hold on;

plot(0:N1,closedloop1.X(1,:),'-','color','blue','LineStyle','-.','LineWidth',3);
hold on;
plot(0:N1,closedloop1.X(2,:),'-','color','yellow','LineStyle',':','LineWidth',3);
hold on;
stairs(0:N1-1,closedloop1.U,'-','color','red','LineStyle','--','LineWidth',3);
hold on;
legend('MPC x_1','MPC x_2','MPC u','closed x_1','closed x_2','closed u');
subplot(2,1,2);
XN1.plot('color','lightblue','alpha',0.6);
hold on;
plot(x0(1),x0(2),'o','MarkerFaceColor','auto');
    plot(data1.X(1,:),data1.X(2,:),'color','k','LineWidth',2)

legend('feasible initial set','initial state','state trajectory');
disp('mpc1 cost=');
disp(data1.cost);

%mpc2
mpc2=MPCController(sys,N2);
[~,feasible,openloop2]=mpc2.evaluate(x0);
[~,feasible,closedloop2]=mpc2.evaluate(x0);

    if feasible==0
        disp('not feasible')
        
    end
loop2 = ClosedLoop(mpc2,sys);
data2=loop2.simulate(x0,N2);
figure(7)
subplot(2,1,1);
stairs(0:N2-1,data2.U,'MarkerSize',12,'LineWidth',1);
hold on;grid on;
plot(data2.X(1,:),'MarkerSize',12,'LineWidth',1);
hold on;
plot(data2.X(2,:),'MarkerSize',12,'LineWidth',1);
xlabel('iteration');
% legend('MPC x_1','MPC x_2','MPC u');
title('MPC2 state trajectories VS prediction');

hold on;
plot(0:N2,closedloop2.X(1,:),'-','color','blue','LineStyle','-.','LineWidth',3);
hold on;
plot(0:N2,closedloop2.X(2,:),'-','color','yellow','LineStyle',':','LineWidth',3);
hold on;
stairs(0:N2-1,closedloop2.U,'-','color','red','LineStyle','--','LineWidth',3);
hold on;
legend('MPC x_1','MPC x_2','MPC u','closed x_1','closed x_2','closed u');
subplot(2,1,2)
XN2.plot('color','lightblue','alpha',0.6);
hold on;
plot(x0(1),x0(2),'o','MarkerFaceColor','auto');
plot(data2.X(1,:),data2.X(2,:),'color','k','LineWidth',2);
% plot(data2.X(end,:),data2.X(end,:),'o','MarkerFaceColor','auto');

legend('feasible initial set','initial state','state trajectory');
disp('mpc2 cost=');
disp(data2.cost);

%mpc3
mpc3=MPCController(sys,N3);
[~,feasible,openloop3]=mpc3.evaluate(x0);
[~,feasible,closedloop3]=mpc3.evaluate(x0);

    if feasible==0
        disp('not feasible')
        
    end
loop3 = ClosedLoop(mpc3,sys);
data3=loop3.simulate(x0,N3);
figure(8)
subplot(2,1,1);
stairs(0:N3-1,data3.U,'MarkerSize',12,'LineWidth',1);
hold on;grid on;
plot(data3.X(1,:),'MarkerSize',12,'LineWidth',1);
hold on;
plot(data3.X(2,:),'MarkerSize',12,'LineWidth',1);
xlabel('iteration');
% legend('MPC x_1','MPC x_2','MPC u');
title('MPC3 state trajectories VS closed loop prediction');

hold on;
plot(0:N3,closedloop3.X(1,:),'-','color','blue','LineStyle','-.','LineWidth',3);
hold on;
plot(0:N3,closedloop3.X(2,:),'-','color','yellow','LineStyle','-','LineWidth',3);
hold on;
stairs(0:N3-1,closedloop3.U,'-','color','red','LineStyle','-','LineWidth',3);
hold on;
legend('MPC x_1','MPC x_2','MPC u','closed x_1','closed x_2','closed u');

subplot(2,1,2);
XN3.plot('color','lightblue','alpha',0.6);
hold on;
plot(data3.X(1,:),data3.X(2,:),'color','k','LineWidth',2)
plot(x0(1),x0(2),'o','MarkerFaceColor','auto');
legend('feasible initial set','initial state','state trajectory');
disp('mpc3 cost=');
disp(data3.cost);






%% Q1c
% find the explicit MPC solution to reach the target set and plot the
% state-space partitions

[P,~,~] = dare(A, B, Q, R);
sys.x.with('terminalPenalty');
sys.x.terminalPenalty = QuadFunction(P);
sys.x.with('terminalSet');
sys.x.terminalSet = Polyhedron('A', eye(2), 'b', 0.01*ones(2,1));


empc3 = mpc3.toExplicit();
eloop3 = ClosedLoop(empc3, sys);
edata3 = eloop3.simulate(x0, N3);
figure(9);
empc3.partition.plot();
xlabel('x1');
ylabel('x2');
title('State-Space Partitions');

%% Q1d
% choose a target set Xf such that persistent feasibility is guaranteed for
% all initial state x0 in C_inf
N_=1;
C_inf = sys.invariantSet();
mpc_=MPCController(sys,N_);
% sys.x.with('terminalSet');
% sys.x.terminalSet = C_inf;
[u, feasible] = mpc_.evaluate(x0);
if ~feasible
    disp('not feasible')
end
Xtarget = sys.reachableSet('X',C_inf,'U',U,'N',1,'direction','forward').intersect(C_inf);
if C_inf.contains(x0)==1
    disp('C_inf contain x0')
end
if Xtarget.contains(x0)==1
    disp('Xtargetset contain x0')
end
figure(10)
X.plot('color','lightblue','alpha',0.5);
hold on;
Xtarget.plot('color','red','alpha',0.6);
hold on;
C_inf.plot('color','yellow','alpha',0.7);
xlabel('x1');ylabel('x2');
legend('X','Xtarget','C\_inf');
Xf_1d = C_inf;%'target set that guarantees persistent feasibility for all initial states inside C_inf';

%% Q2a
% find the matrices for the discretized system

% Define the system parameters
R = 20; % Resistance of armature
kT = 10; % Motor constant
p = 20; % Gear ratio
kJ = 1280.2; % Torsional rigidity
JM = 0.5; % Motor inertia
BM = 0.1; % Motor viscous friction coefficient
JL = 50*JM; % Nominal load inertia
BL = 25; % Load viscous friction coefficient

% Define the continuous-time state-space matrices
A_c = [0 1 0 0;
       -kJ/JL -BL/JL kJ/p/JL 0;
       0 0 0 1;
       kJ/p/JM 0 -kJ/p^2/JM -(BM+kT^2/R)/JM];
B_c = [0;
       0;
       0;
       kT/R/JM];
C_c = [kJ 0 -kJ/p 0];

% Sampling interval
h = 0.1;

% Discretize the system
sys_d = c2d(ss(A_c, B_c, C_c, 0), h, 'zoh');

% Extract the discrete-time matrices
A_d = sys_d.A;
B_d = sys_d.B;
C_d = sys_d.C;

% Display the matrices
disp('Discrete-time A matrix:');
disp(A_d);
disp('Discrete-time B matrix:');
disp(B_d);
disp('Discrete-time C matrix:');
disp(C_d);


A_2a =A_d;% 'A matrix for the discretized system';
B_2a = B_d;%'B matrix for the discretized system';
C_2a = C_d;%'C matrix for the discretized system';

%% Q2b
% design a minimum time (i.e. minimum N) controller that brings the system
% to standstill and plot the predicted states and the output
% Minimum-time control problem setup
x0 = [0; 2.5; 0; 75]; % Initial state
xf = [0; 0; 0; 0]; % Target state
u_max = 200; % Maximum input voltage
T_max = 150; % Maximum torsional torque

% Initializations
x = x0;
time_steps = 0;
u = [];


% Time settings
h = 0.1;  % Sampling interval
t_final = 20;  % Assume a final time, adjust based on system response
t = 0:h:t_final;  % Time vector

% Initialize variables
x = zeros(size(A_d, 1), length(t));  % State trajectory
x(:,1) = x0;  % Set initial state
u = zeros(1, length(t));  % Control input trajectory

% Control loop
for k = 1:(length(t) - 1)
    % Calculate the output torque T for the current state
    T = C_d * x(:,k);
    
    % Check if the output torque is within the limits
    if abs(T) > T_max
        % If the torque exceeds the maximum, saturate it
        T = sign(T) * T_max;
    end
    
    % Decide on control action
    if x(1,k) > 0 || x(2,k) > 0  % Assuming we want to reach zero position and velocity
        u(k) = -u_max;  % Apply negative maximum voltage to decelerate
    else
        u(k) = u_max;   % Apply positive maximum voltage to accelerate
    end
    
    % Ensure control input does not exceed maximum voltage
    if abs(u(k)) > u_max
        u(k) = sign(u(k)) * u_max;
    end
    
    % Update the state using the system dynamics
    x(:,k+1) = A_d * x(:,k) + B_d * u(k);
    
    % If the system is close to the desired state, break the loop
    if all(abs(x(:,k+1)) < 1e-3)  % Threshold of 1e-3 for reaching the target state
        break;
    end
end

% Plot the results
figure;
subplot(2,1,1);
plot(t, x,'LineWidth', 2);
title('State Trajectories');
xlabel('Time (s)');
ylabel('States');
legend('theta_L', 'dtheta_L/dt', 'theta_M', 'dtheta_M/dt');

subplot(2,1,2);
stairs(t, u, 'b', 'LineWidth', 2); % Plot the control input
hold on; % Hold the plot for superimposing the next signal
stairs(t, T_max*ones(size(t)), 'r--', 'LineWidth', 2); % Plot the torque limit as a red dashed line
hold off; % Release the plot hold

title('Control Input and Torque Limit');
xlabel('Time (s)');
ylabel('Value');
legend('Control input u', 'Torque limit T_{max}');
grid on; % Add grid for better readability


% Assume x is your state trajectory matrix where
% x(2,:) corresponds to theta_L_dot and x(4,:) corresponds to theta_M_dot

% Define a convergence threshold (this could be a very small number, like 1e-3)
threshold = 1e-3;

% Find the first index where both velocities are below the threshold
convergence_idx = find(abs(x(2,:)) < threshold & abs(x(4,:)) < threshold, 1, 'first');

% Check if a convergence index was found
if isempty(convergence_idx)
    disp('System did not converge within the simulated time.');
else
    % Calculate the minimum time to convergence
    minimum_time = convergence_idx * h;  % h is your sampling time

    % Plotting the state trajectories
    figure;
    plot(t, x(2,:), t, x(4,:));  % Plot theta_L_dot and theta_M_dot
    hold on;

    % Mark the convergence point on the plot
    plot(t(convergence_idx), x(2,convergence_idx), 'ro');  % Mark theta_L_dot convergence
    plot(t(convergence_idx), x(4,convergence_idx), 'ro');  % Mark theta_M_dot convergence
    title('State Trajectories with Convergence Point');
    xlabel('Time (s)');
    ylabel('State Values');
    legend('\theta_{L_{dot}}', '\theta_{M_{dot}}', 'Convergence Point');
    grid on;

    % Annotate the minimum time on the plot
    xline(t(convergence_idx), '--', sprintf('Minimum Time: %.2f s', minimum_time));
    hold off;
end








%% Q2c
% design a minimum time controller that brings the system to standstill,
% while having a constrained torsional torque, and plot the predicted states
% and the output

N_2c = %'minimum horizon length that allows the system to reach standstill while respecting the constraint on T';

%% Q2d
% design a minimum time controller that brings the system close to
% standstill and plot the predicted states and output

% then, perform a closed loop simulation with a controller having the N you
% just found and plot the states and outputs for 2s

N_2d = %'minimum horizon length that allows the system to reach almost standstill';