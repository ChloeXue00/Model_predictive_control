clc
close all
clear

%% Q1a
% represent the polyhedron in V- and H-representation and plot the two
% versions separately


%H-representation:
A=[0,1;
   -1,0;
   -1,-1; 
   1,1];
b=[0;
   0;
   1;
   1];

P1 = Polyhedron(A, b);

disp('Polyhedron P1:');
axis on; grid on; hold on; 
subplot(1,2,1);
title('Polyhedron P1');
plot(P1,'color','lightgreen','alpha',0.4); axis([-2 2 -2 2]);
xlabel('$x_1$'); ylabel('$x_2$'); 
hold on; 


 %V representation
V1 = [0 0]; V2 = [1 0]; V3 = [0 -1]; V4 = [1.41421 -2.41421];V5=[2.41421,-1.41421];
V=[V1;V2;V3;V4;V5];
V1 = Polyhedron(V);
 % This plots the vertices as red circles
disp('Polyhedron V1:');
axis on; grid on; hold on; 

subplot(1,2,1);
plot(V1,'color','lightgreen','alpha',0.4);
title('Polyhedron P1');
axis([-2 2 -2 2]);
xlabel('$x_1$'); ylabel('$x_2$');  hold off; 



subplot(1,2,2);
plot(V1,'color','lightgreen','alpha',0.4);
title('Polyhedron V1');
axis([-2 2 -2 2]);
xlabel('$x_1$'); ylabel('$x_2$'); hold off; 


P_H_1a = P1; %'H-representation of the polyhedron';
P_V_1a = V1; %'V-representation of the polyhedron';

%% Q1b
% perform some sets operations
A1=[0,1;
    1,0;
    0,-1;
    -1,0];
b1 = [2;2;2;2];
P1b = Polyhedron(A1, b1);

A2=[-1,-1;
    1,1;
    1,-1;
    -1,1];
b2=[1;1;1;1];
Q1b=Polyhedron(A2, b2);



P_1b = P1b;%'Polyhedron P';

Q_1b =Q1b;% 'Polyhedron Q';

Mink_sum_1b = P1b + Q1b; %'P + Q';

Pontr_diff_1b = P1b-Q1b; %'P - Q';

PminusQplusQ_1b = (P1b-Q1b)+Q1b;%   '(P - Q) + Q';

PplusQminusQ_1b =(P1b + Q1b)-Q1b; % '(P + Q) - Q';

QminusPplusP_1b = (P1b-Q1b)+P1b;  % '(Q - P) + P';

QplusPminusP_1b = (P1b + Q1b)-P1b; % '(Q + P) - P';



figure;
subplot(4,2,1);
plot(P1b, 'color', 'lightblue','alpha',0.4);
xlabel('$x_1$'); ylabel('$x_2$'); 
% title('Polyhedron P');
title('P');

subplot(4,2,2);
plot(Q1b, 'color', 'lightblue','alpha',0.4);
xlabel('$x_1$'); ylabel('$x_2$'); 
% title('Polyhedron Q');
title('Q');


subplot(4,2,3);
plot(Mink_sum_1b,'color', 'lightyellow','alpha',0.4);
xlabel('$x_1$'); ylabel('$x_2$'); 
% title('Mink\_sum\_1b');
title('P + Q');


subplot(4,2,4);
plot(Pontr_diff_1b,'color', 'lightgreen','alpha',0.4);
xlabel('$x_1$'); ylabel('$x_2$'); 
title('P - Q');

subplot(4,2,5);
plot(PminusQplusQ_1b,'color', 'lightyellow','alpha',0.4);
xlabel('$x_1$'); ylabel('$x_2$'); 
title('(P- Q) + Q');


subplot(4,2,6);
plot(PplusQminusQ_1b,'color', 'lightgreen','alpha',0.4);
xlabel('$x_1$'); ylabel('$x_2$'); 
title('(P+Q) - Q');


subplot(4,2,7);
plot(QminusPplusP_1b,'color', 'lightyellow','alpha',0.4);
xlabel('$x_1$'); ylabel('$x_2$'); 
title('(Q-P)+P');



subplot(4,2,8);
plot(QplusPminusP_1b,'color', 'lightgreen','alpha',0.4);
xlabel('$x_1$'); ylabel('$x_2$'); 
title('(Q+P)-P');



%% Q2a
% show that the set S is positively invariant for the system and plot it
A_d= [0.8,0.4;
    -0.4,0.8];  %x(k+1) = A*x(k);
%S={Ain x <=bin};
A3=[1,0;
    0,1;
    -1,0;
    0,-1;
    1,1;
    1,-1;
    -1,1;
    -1,-1];
b3 = [1;1;1;1;1.5;1.5;1.5;1.5];
% %find x for inequality
% options = optimoptions('linprog', 'Display', 'none');
% 
% % Since we only want to find a feasible point, the objective function can be zero
% f = [0; 0];

% % Use linprog to find a point x that satisfies A_in*x ≤ b_in
% [x, fval, exitflag] = linprog(f, A_in, b_in, [], [], [], [], options);
% % Check if a feasible point was found
% if exitflag == 1
%     % A feasible point x was found
%     disp('A point that satisfies A_in*x ≤ b_in is:');
%     disp(x);
% else
%     % No feasible point was found
%     disp('No feasible point was found.');
% end
% Check if the polyhedron is bounded (in which case it has vertices)

% P3 = Polyhedron(A3, b3);
% 

S = Polyhedron('A', A3, 'b', b3);

AS = S.affineMap(A_d);

if AS <= S
    disp('The set S is positively invariant.');
else
    disp('The set S is not positively invariant.');
end



S.plot('color', 'lightblue','alpha',0.4);
hold on;
AS.plot('color', 'lightblue','alpha',0.7);
legend('off');
legend('S', 'AS');
% legend('S','AS');

S_2a = S; % 'Set S';

%% Q2b
% calculate the one-step reachable set from S and plot it together with S

A_d= [0.8,0.4;
    -0.4,0.8];  %x(k+1) = A*x(k);
B_d = [0;1];


U = Polyhedron('lb', -1, 'ub', 1);
lb=-1;
ub=1;
reachablePoints = [];

%affine mapping
A_mapped_X = A_d * S;

% Compute the Minkowski sum of the resulting polytopes (A_mapped_X and B*U)
B_mapped_U = B_d * U;
reach_set = plus(A_mapped_X, B_mapped_U);

% Plot the original set S and the reachable set with transparency
figure;
plot(S,'color', 'blue', 'alpha', 0.5);
hold on;
plot(reach_set, 'color', 'lightblue','alpha', 0.5);
xlabel('$x_1$');
ylabel('$x_2$');
title('Original Set S and One-step Reachable Set');
legend('Set S', 'Reach Set');



 

reach_set_S_2b =reach_set;% 'One-step reachable set from S';

%% Q2c
% calculate the one-step precursor set of S and plot it together with S

Ain = [1,0;
    0,1;
    -1,0;
    0,-1;
    1,1;
    1,-1;
    -1,1;
    -1,-1];
bin=[1;
    1;
    1;
    1;
    1.5;
    1.5;
    1.5;
    1.5];

bin2=[1;1];
Ain2=[1 ;-1]; 
A_aug = [Ain*A_d, Ain*B_d; zeros(2,2), Ain2];
b_aug = [bin;bin2];

S_pre = Polyhedron('A',A_aug,'b',b_aug);

% projecting the preS set onto the x-space

Pre_prj=projection(S_pre,(1:2));


figure;

S_pre.plot('color', 'lightblue', 'alpha', 0.5);
hold on;
Pre_prj.plot('color', 'yellow', 'alpha', 0.5)
legend({'PreSet S', 'PreSet projection'});
hold off;




pre_set_S_2c = Pre_prj;%'One-step precursor set of S'; % project tau into 1st and 2nd dimensions






%% Q3a ,parameters setting
% find shortest prediction horizon such that the RH controller is feasible


% Define the system matrices
A = [0.9 0.4; -0.4 0.9];
B = [0; 1];


Q=eye(2);
R=1;

x0 = [2;0];

sys = LTISystem('A', A, 'B', B);     
sys.x.min = [-3; -3];
sys.x.max = [ 3;  3];
sys.u.min = -0.1;
sys.u.max =  0.1;

X = Polyhedron('lb', sys.x.min, 'ub', sys.x.max);
U = Polyhedron('lb', sys.u.min, 'ub', sys.u.max);

sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);

% terminal cost
Pf=zeros(2);
sys.x.with('terminalPenalty');
sys.x.terminalPenalty = QuadFunction(Pf);

% terminal set
Xf   = Polyhedron([0 0]);
sys.x.with('terminalSet');
sys.x.terminalSet = Xf;

% Q3a, loop to find optimal N

% Initialize horizon
N = 1;
max_horizon = 30; % Define a reasonable upper limit for N

% Initialize control system
mpc = [];

% Loop to find the shortest feasible horizon
while N <= max_horizon
    % Attempt to construct the MPC controller
    try
        mpc = MPCController(sys, N);
        sys.x.terminalSet = Xf;

        % Check feasibility for the initial condition
        x0 = [2; 0];
        [u, feasible, openloop] = mpc.evaluate(x0);
        
        % If a feasible solution is found, exit the loop
        if feasible
            fprintf('Feasible solution found for N = %d\n', N);
            break;
        end
    catch ME
        % If an error occurs (e.g., the problem is infeasible), continue to the next N
        disp(ME.message);
    end
    
    % Increment the horizon
    N = N + 1;
end



C_inf1  = sys.invariantSet(); % largest control invariant set


N_3a = 26;%'Shortest prediction horizon that makes the RH controller feasible';

%% Q3b
% find out if the RH controller is persistently feasible all the way to the
% origin

% Set N = 2 for the prediction horizon
N2 = 2;

C_inf2  = sys.invariantSet(); % largest control invariant set


% Define the terminal set as the origin

Xf2 = C_inf2;

% Xf1   = Polyhedron([0 0]);
% sys.x.with('terminalSet');
sys.x.terminalSet = Xf2;

mpc_2 = MPCController(sys, N2);
[~, feasible, ~] = mpc_2.evaluate(x0);
if feasible
        fprintf('The RH controller is persistently feasible for all tested initial states.\n');
else
    fprintf('The RH controller is not persistently feasible. Tested %d initial states, %d were feasible.\n', size(initialStates, 2), feasibleCount);
end


%The RH controller is persistently feasible for all tested initial states.


maxContrInvSet_3b = C_inf2;%'Maximal control invariant set for the system';

%% Q3c
% plot the set of feasible initial states for the two controllers
% previously designed and discuss about the size of the problems

% Calculate the feasible set for the controller setup
XN_first = Xf;
N1=26;
for i = 1:N1
    XN_first = sys.reachableSet('X', XN_first, 'U', U, 'direction', 'backward');
    XN_first = XN_first.intersect(X).minHRep();
end


XN_second = Xf2;
N2=2;
for i = 1:N2
    XN_second = sys.reachableSet('X', XN_second, 'U', U, 'direction', 'backward');
    XN_second = XN_second.intersect(X).minHRep();
end

% Plot the feasible set for the first controller setup
figure; 
XN_first.plot('color', 'blue', 'alpha', 0.5);
hold on;
XN_second.plot('color', 'lightblue', 'alpha', 0.5);
title('Feasible Initial States for first and Second RH Controller');
xlabel('State 1'); ylabel('State 2');
legend('feasible set for first controller: X_f = 0 and shortest N=26', ...
    'feasible set for second controller: X_f = C_{inf},N=2');










X0_1_3c =XN_first;% 'Set of feasible initial states for controller designed in Q3a';
X0_2_3c = XN_second;%'Set of feasible initial states for controller designed in Q3b';
