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
k5 = zeros(1,n+1);
k6 = zeros(1,n+1);
y(1) = ya;

for i = 1:1:n+1
    x(i) = a + (i-1)*h;
end

for i = 1:1:n
   k1(i) = f(x(i),y(i));
   k2(i) = f(x(i) + h/4, y(i) + (k1(i)*h)/4);
   k3(i) = f(x(i) + h/4, y(i) + (k1(i)*h)/8 + (k2(i)*h)/8);
   k4(i) = f(x(i) + h/2, y(i) - (k2(i)*h)/2 + k3(i)*h);
   k5(i) = f(x(i) + (3*h)/4, y(i) + (3*k1(i)*h)/16 + (9*k4(i)*h)/16);
   k6(i) = f(x(i) +h, y(i) - (3/7)*(k1(i)*h) + (2/7)*(k2(i)*h) + (12/7)*(k3(i)*h) - (12/7)*(k4(i)*h) + (8/7)*(k5(i)*h));
   y(i+1) = y(i) + (h/90)*(7*k1(i) + 32*k3(i) + 12*k4(i) + 32*k5(i) + 7*k6(i));
end

disp("k1 k2 k3 k4 k5 k6 xi yi - columns respectively")
k = [transpose(k1) transpose(k2) transpose(k3) transpose(k4) transpose(k5) transpose(k6)];
disp(k);
result = [transpose(x) transpose(y)];
disp("x y - columns respectively")
disp(result);
% Unccoment the following line to plot graph
% plot(x,y);

function fx = f(x,y)
    % Enter your function here. If the function doesnt depend on y then just
    % add an extra term '0*y' at the end to avoid getting an error
    fx = 4*exp(0.8*x)-0.5*y;
end