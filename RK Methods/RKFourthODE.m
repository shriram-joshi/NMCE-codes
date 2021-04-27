% By Shriram Joshi
clc;
clearvars;
% Change to 'short' or 'short e' for truncated decimals
format long;

% Input range, step size and y(a) i.e your initial guess of y
a = 0;
b = 0.5;
h = 0.5;
ya = 2;

% Calculations
n = (b-a)/h;
x = zeros(1,n+1);
y = zeros(1,n+1);
k1 = zeros(1,n+1);
k2 = zeros(1,n+1);
k3 = zeros(1,n+1);
k4 = zeros(1,n+1);
y(1) = ya;

for i = 1:1:n+1
    x(i) = a + (i-1)*h;
end

for i = 1:1:n
   k1(i) = f(x(i),y(i));
   k2(i) = f(x(i) + h/2, y(i) + (k1(i)*h)/2);
   k3(i) = f(x(i) + h/2, y(i) + (k2(i)*h)/2);
   k4(i) = f(x(i) + h, y(i) + k3(i)*h);
   y(i+1) = y(i) + (h/6)*(k1(i) + 2*k2(i) + 2*k3(i) + k4(i));
end

disp("k1 k2 k3 k4 xi yi - columns respectively")
result = [transpose(k1) transpose(k2) transpose(k3) transpose(k4) transpose(x) transpose(y)];
disp(result);
% Unccoment the following line to plot graph
% plot(x,y);

function fx = f(x,y)
    % Enter your function here. If the function doesnt depend on y then just
    % add an extra term '0*y' at the end to avoid getting an error
    fx = 4*exp(0.8*x)-0.5*y;
end