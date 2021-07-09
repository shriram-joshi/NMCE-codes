% The Crank - Nicholson Method for Transient Heat Conduction in an
% insulated rod. This code is using non-dimensional method and the
% delta_t_bar increases by 10% for every iteration.
% By Shriram Joshi

clc;
clearvars;
format short;

% Inputs
% k = k'/rho*C where k' is thermal conductivity
k = 0.835;
rho = 2.7;
C = 0.2174;
L = 10;
Dx = 2;
Dt = 1;
% Number of iterations
iter= 10;

Dx_bar = Dx/L;
Dt_bar = k*Dt/(rho*C*L^2);

% Enter inital values of internal nodes in this array. The first element 
% can be considered T_0 and last element as T_L. Ensure the length matches
% with the number of nodes formed by the entered delta x 
T_init = [100, 0, 0, 0, 0, 50];

% Number of internal nodes
m = L/Dx - 1;

if(length(T_init) ~= m+2)
    disp("Error: The number of nodes")
    disp("(in this case " + length(T_init) + ") in your T_init matrix does not match")
    disp("with the number of nodes (in this case " + (m+2) + ")")
    disp("that are formed if we consider delta_x as " + Dx + " and L as " + L)
    return
end

A = zeros(m,m);
B = zeros(m,1);

% We set dimension of matrix u to iter+1 because we are including the time instance t = 0 
%  and m+2 because there are m internal and 2 end nodes
u = zeros(iter+1,m+2);


for i = 1:length(T_init)
    u(1,i) = (T_init(i) - T_init(1)) / (T_init(m+2) - T_init(1));
end

for i = 2:iter+1
    % Recalculating lambda after increasing Dt_bar by 10%
    lambda = rho*C*Dt_bar/Dx_bar^2;
    
    % Initializing the A matrix in [A][u] = [B]
    for j = 1:m
        A(j,j) = 2*(1+lambda);
        if(j~=m)
            A(j+1,j) = -lambda;
            A(j,j+1) = -lambda;        
        end
    end
    
    
    % Initializing the B matrix in [A][u] = [B] 
    % for every iteration (t_bar = n*Dt_bar, n = 1, 2, 3,... ,iter)
    B(1) = 2*lambda*u(i-1,1) + 2*(1-lambda)*u(i-1,2) + lambda*u(i-1,3);
    for j = 2:m-1
        B(j) = lambda*u(i-1,j) + 2*(1-lambda)*u(i-1,j+1) + lambda*u(i-1,j+2);
    end
    B(m) = lambda*u(i-1,m) + 2*(1-lambda)*u(i-1,m+1) + 2*lambda*u(i-1,m+2);
    
    % Solve for internal node values
    u_internal = linsolve(A,B);
    
    % Varies the values of end nodes as per some formula
    % Edit if values of end nodes change with time
    u(i,1) = u(i-1,1);
    u(i,m+2) = u(i-1,m+2);
    
    % Enters the values of internal nodes at time t_bar into the u matrix
    for j = 2:m+1
        u(i,j) = u_internal(j-1); 
    end
    
    % Increasing Dt_bar by 10% at every iteration for faster calculation
    Dt_bar = Dt_bar + 0.1*Dt_bar;
end


% Initializing x_bar values according to Dx_bar for ploting
x_bar = zeros(m+2,1);
for i= 2:m+2
    x_bar(i) = x_bar(i-1) + Dx_bar;
end

% Plots u vx x_bar for different t_bar
% figure("name","u vs x_bar")
plot(x_bar,u);
grid;

% % Initializing x and T values for ploting
% x = x_bar.*L;
% T = zeros(iter+1,m+2);
% T(1,:) = T_init;
% T = u.*(T(1,m+2)-T(1,1)) + T(1,1);
% 
% % Plots T vx x for different t
% figure("name", "T vs x (via Non-Dimensional)")
% plot(x,T);
% grid;