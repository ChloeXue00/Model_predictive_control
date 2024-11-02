clc
close all
clear all

%% Q3c
% solve the LP using the formulation shown in (5)
% Define the matrix A and vector b as given in the problem
A = [0.4889, 0.2939;
     1.0347, -0.7873;
     0.7269, 0.8884;
     -0.3034, -1.1471];
b = [-1.0689;
     -0.8095;
     -2.9443;
     1.4384];

% The problem to solve is min ||Ax - b||_infinity
% This can be written as a linear program where we minimize epsilon subject to
% -epsilon <= Ax - b <= epsilon

% Number of constraints
n = size(A,2);

% Construct the matrices for the constraints

F = [A, -ones(length(b), 1); 
     -A, -ones(length(b), 1)];

g = [b; -b];


% Objective function: we want to minimize epsilon, so the cost vector c is 
c = [zeros(n, 1);1]; % Coefficients for x are 0, coefficient for epsilon is 1

% Lower bounds (no lower bound for epsilon and x)
lb = [-inf(n, 1); -inf];
% No upper bound for all variables
ub = [];
% Use linprog to solve the linear program
options = optimoptions('linprog','Algorithm','dual-simplex');
[x, fval] = linprog(c, F, g, [], [], lb, [], options);

% The optimal value of epsilon is the last entry in the solution vector x
epsilon_optimal = x(3);
% The solution vector x is the first entries in the solution vector x
x_optimal = x(1:2);
disp('optimal epsilon=')
disp(x(3));
disp('optimal x=')
disp(x(1:2));
z_3c = x(3);%'Solution to the linear problem';

%% Q3e
A = [0.4889   0.2939;
     1.0347  -0.7873;
     0.7269   0.8884;
     -0.3034  -1.1471];

b = [-1.0689;
     -0.8095;
     -2.9443;
      1.4384];




% Dual problem
F = [A,-ones(4,1);
    -A,-ones(4,1)];
Aeq= F';
beq= -c;
f = [b; -b];  
lb = zeros(length(f),1);
ub = [];
y=linprog(f,[],[],Aeq,beq,lb,ub);

mu_3e = y;%'Solution to the dual problem';

%% Q3f
% Solve the primal problem by using the solution of the dual problem
% obtained in Q3e


c = [0; 0; 1];
F = [0.4889, 0.2939, -1.0000;
     1.0347, -0.7873, -1.0000;
     0.7269, 0.8884, -1.0000;
     -0.3034, -1.1471, -1.0000;
     -0.4889, -0.2939, -1.0000;
     -1.0347, 0.7873, -1.0000;
     -0.7269, -0.8884, -1.0000;
     0.3034, 1.1471, -1.0000];
g = [-1.0689;
     -0.8095;
     -2.9443;
     1.4384;
     1.0689;
     0.8095;
     2.9443;
     -1.4384];

% Define the dual variables from the dual solution
y = [0; 0; 0.4095; 0.4284; 0; 0.1621; 0; 0];

% Solve the primal problem using linprog
options = optimoptions('linprog','Algorithm','dual-simplex');
[z, fval, exitflag, output] = linprog(c, F, g, [], [], [], [], options);
% Display the optimal solution to the primal problem
disp('Optimal z:');
disp(z);


z_3f =z; % 'Solution to the primal problem obtained through the dual one';

%% Q4a
% solve the QP


% Define the H matrix for the quadratic terms in the objective function
H = 2*eye(4); 

% Define the f vector for the linear terms in the objective function
f = zeros(4,1);

% No linear inequalities or equalities
A = 0.4;
b = 1;


% Define the bounds on the variables x1, x2, u0, u1
lb = [2.5; -0.5; -2; -2];
ub = [5; 0.5; 2; 2];

% Initial condition
x0 = 1.5;

% The state transition matrix for x1 and x2 based on the system dynamics (assuming a two-step process)
% x1 = 0.4*x0 + u0
% x2 = 0.4*x1 + u1

%x=[x1,x2,u0,u1]
% State transition constraints as equality constraints
Aeq = [1 , 0 , -b , 0;
      -A ,1  ,0 ,-b];
beq = [A*x0; 0];

% Solve the quadratic programming problem
options = optimoptions('quadprog','Display','off');
%[x, fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
%[x, fval,~, output, lambda] = quadprog(H, f, A, b, Aeq, beq, lb, ub, x0, options);
[x, fval,~, output, lambda] = quadprog(H, f, [], [], Aeq, beq, lb, ub, x0, options);

% Display the solution
disp('Solution:');
disp(['x1: ', num2str(x(1))]);
disp(['x2: ', num2str(x(2))]);
disp(['u0: ', num2str(x(3))]);
disp(['u1: ', num2str(x(4))]);
disp(['Objective function value: ', num2str(fval)]);
disp(lambda);


x_4a =[x(1);x(2)];% 'Solution to the QP (x1 and x2)';
u_4a =[x(3);x(4)];% 'Solution to the QP (u0 and u1)';


%  lambda is the structure returned by quadprog

% Check for active lower bound constraints
active_lower_bounds = find(lambda.lower > 1e-6);  % 1e-6 is a small tolerance

% Check for active upper bound constraints
active_upper_bounds = find(lambda.upper > 1e-6);  % 1e-6 is a small tolerance

% Check for inequality upper bound constraints
active_ineqlin_bounds = find(lambda.ineqlin> 1e-6);

% Display active constraints
disp('Active lower bound constraints at indices:');
disp(active_lower_bounds);

disp('Active upper bound constraints at indices:');
disp(active_upper_bounds);

disp('Active inequality bound constraints at indices:');
disp(active_ineqlin_bounds);

% Since lambda_eqlin is negative , meaning the optimal solution is 
%not reaching the equality constraint, which is reasonable in this problem.
% active_eqlin_bounds = find(lambda.eqlin> 1e-6);
% disp('Active equality bound constraints at indices:');
% disp(active_eqlin_bounds);  