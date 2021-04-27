% By Shriram Joshi
clc;
clearvars;
format long;

% Input range, step size and y(a) i.e your initial guess of y
a = 0;
b = 4;
h = 0.5;
ya = 1;

% Calculations
n = (b-a)/h;
x = zeros(1,n+1);
y = zeros(1,n+1);
y(1) = ya;

for i = 1:1:n+1
    x(i) = a + (i-1)*h;
end

for i = 1:1:n
    y(i+1) = y(i) + f(x(i),y(i))*h;
end

disp("x and y")
result = [transpose(x) transpose(y)];
disp(result);
% Unccoment the following line to plot graph
% plot(x,y);

function fx = f(x,y)
    % Enter your function here. If the function doesnt depend on y then just
    % add an extra term '0*y' at the end to avoid getting an error
    fx = -2*x.^3 + 12*x.^2 - 20*x + 8.5 + 0*y;
end