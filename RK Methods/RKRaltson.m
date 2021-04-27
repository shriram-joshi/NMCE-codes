% By Shriram Joshi
clc;
clearvars;
% Change to 'short' or 'short e' for truncated decimals
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
k1 = zeros(1,n+1);
k2 = zeros(1,n+1);
y(1) = ya;

for i = 1:1:n+1
    x(i) = a + (i-1)*h;
end

for i = 1:1:n
   k1(i) = f(x(i),y(i));
   k2(i) = f(x(i)+(3*h)/4, y(i) + (3*k1(i)*h)/4);
   y(i+1) = y(i) + (k1(i)/3 + (2*k2(i))/3)*h;
end

disp("k1 k2 x y - columns respectively")
result = [transpose(k1) transpose(k2) transpose(x) transpose(y)];
disp(result);
% Unccoment the following line to plot graph
% plot(x,y);

function fx = f(x,y)
    % Enter your function here. If the function doesnt depend on y then just
    % add an extra term '0*y' at the end to avoid getting an error
    fx = -2*x.^3 + 12*x.^2 - 20*x + 8.5 + 0*y;
end