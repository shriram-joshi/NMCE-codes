% The Crank - Nicholson Method for Transient Heat Conduction in an
% insulated rod.
% By Shriram Joshi

clc;
clearvars;
format short;

% Inputs
% k = k'/rho*C where k' is thermal conductivity
k = 0.835;
L = 10;
Dx = 2;
Dt = 1;
% Number of iterations
iter= 10;

% Enter inital values of internal and external nodes in this array, Ensure it matches
% with the number of nodes formed by the entered delta x 
T_init = [100, 0, 0, 0, 0, 50];

lambda = k*Dt/Dx^2;
% Number of internal nodes
m = L/Dx - 1;

if(length(T_init) ~= m+2)
    disp("Error: The number of internal nodes")
    disp("(in this case " + length(T_init)-2 + ") in your T_init matrix does not match")
    disp("with the number of internal nodes (in this case " + m + ")")
    disp("that are formed if we consider delta_x as " + Dx + " and L as " + L)
    return
end

A = zeros(m,m);
B = zeros(m,1);

% We set dimension of matrix T to iter+1 because we are including the time instance t = 0 
%  and m+2 because there are m internal and 2 end nodes
T = zeros(iter+1,m+2); 

T(1,:) = T_init;

% Initializing the A matrix in [A][T] = [B]
for i = 1:m
    A(i,i) = 2*(1+lambda);
    if(i~=m)
        A(i+1,i) = -lambda;
        A(i,i+1) = -lambda;        
    end
end


for i = 2:iter+1
    % Initializing the B matrix in [A][T] = [B] 
    % for every iteration (t = n*Dt, n = 1, 2, 3,... ,iter)
    B(1) = 2*lambda*T(i-1,1) + 2*(1-lambda)*T(i-1,2) + lambda*T(i-1,3);
    for n = 2:m-1
        B(n) = lambda*T(i-1,n) + 2*(1-lambda)*T(i-1,n+1) + lambda*T(i-1,n+2);
    end
    B(m) = lambda*T(i-1,m) + 2*(1-lambda)*T(i-1,m+1) + 2*lambda*T(i-1,m+2);
    
    % Solve for internal node values
    T_internal = linsolve(A,B);
    
    % Varies the values of end nodes as per some formula
    % Edit if values of end nodes change with time
    T(i,1) = T(i-1,1);
    T(i,m+2) = T(i-1,m+2);
    
    % Enters the values of internal nodes at time t into the T matrix
    for n = 2:m+1
        T(i,n) = T_internal(n-1); 
    end
end

% Initializing x values according to Dx for ploting
x = zeros(m+2,1);
for i= 2:m+2
    x(i) = x(i-1) + Dx;
end

% Plots T vx x for different t
figure("name","T vs x (Dimensional)");
plot(x,T);
grid;