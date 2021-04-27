% By Shriram Joshi
clc;
clearvars;
% Change to 'short' or 'short e' for truncated decimals
format short;

% Input range, step size and y(a) i.e your initial guess of y
a = 0;
b = 2;
h = 0.5;
y1 = 4;
y2 = 6;

% Calculations
n = (b-a)/h;
x = zeros(1,n+1);
y = zeros(2,n+1);
k1 = zeros(2,n+1);
k2 = zeros(2,n+1);
k3 = zeros(2,n+1);
k4 = zeros(2,n+1);
y(1,1) = y1;
y(2,1) = y2;

for i = 1:1:n+1
    x(i) = a + (i-1)*h;
end

for i = 1:1:n
   [k1(1,i), k1(2,i)] = f(x(i),y(1,i), y(2,i));
   [k2(1,i), k2(2,i)] = f(x(i) + h/2, y(1,i) + (k1(1,i)*h)/2, y(2,i) + (k1(2,i)*h)/2);
   [k3(1,i), k3(2,i)] = f(x(i) + h/2, y(1,i) + (k2(1,i)*h)/2, y(2,i) + (k2(2,i)*h)/2);
   [k4(1,i), k4(2,i)] = f(x(i) + h, y(1,i) + k3(1,i)*h, y(2,i) + k3(2,i)*h);
   y(1,i+1)= y(1,i) + (h/6)*(k1(1,i) + 2*k2(1,i) + 2*k3(1,i) + k4(1,i));
   y(2,i+1)= y(2,i) + (h/6)*(k1(2,i) + 2*k2(2,i) + 2*k3(2,i) + k4(2,i));
end

disp("k11 k12 k21 k22 k31 k32 k41 k42 - columns respectively")
k = [transpose(k1) transpose(k2) transpose(k3) transpose(k4)];
disp(k);
result = [ transpose(x) transpose(y)];
disp("x y1 y2 - columns respectively")
disp(result);

function [fx1, fx2] = f(x,y1,y2)
    % Enter your functions here. If the function doesnt depend on x then just
    % add an extra term '0*x' at the end to avoid getting an error
    fx1 = -0.5*y1 + 0*x;
    fx2 = 4 - 0.3*y2 - 0.1*y1;
end