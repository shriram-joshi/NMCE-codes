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

% Instantiation and Calculations
n = (b-a)/h;
x = zeros(1,n);
y = zeros(2,n);
k1 = zeros(1,2);
k2 = zeros(1,2);
k3 = zeros(1,2);
k4 = zeros(1,2);
y(1,1) = y1;
y(2,1) = y2;

for i = 1:1:n+1
    x(i) = a + (i-1)*h;
end

for i = 1:1:n
   [k1(1), k1(2)] = f(x(i),y(1,i), y(2,i));
   [k2(1), k2(2)] = f(x(i) + h/2, y(1,i) + (k1(1)*h)/2, y(2,i) + (k1(2)*h)/2);
   [k3(1), k3(2)] = f(x(i) + h/2, y(1,i) + (k2(1)*h)/2, y(2,i) + (k2(2)*h)/2);
   [k4(1), k4(2)] = f(x(i) + h, y(1,i) + k3(1)*h, y(2,i) + k3(2)*h);
   y(1,i+1)= y(1,i) + (h/6)*(k1(1) + 2*k2(1) + 2*k3(1) + k4(1));
   y(2,i+1)= y(2,i) + (h/6)*(k1(2) + 2*k2(2) + 2*k3(2) + k4(2));
   
   disp("For iteration " + i + ":");
   for c = 1:1:2
      disp("k1"+c+"="+k1(c));
      disp("k2"+c+"="+k2(c));
      disp("k3"+c+"="+k3(c));
      disp("k4"+c+"="+k4(c));
   end
   
   disp("y1("+x(i+1)+")= " + y(1,i+1));
   disp("y2("+x(i+1)+")= " + y(2,i+1));
   disp("------------------");
end

result = [ transpose(x) transpose(y)];
disp("Final result:-");
disp("x y1 y2 - columns respectively")
disp(result);

function [fx1, fx2] = f(x,y1,y2)
    % Enter your functions here. If the function doesnt depend on x then just
    % add an extra term '0*x' at the end to avoid getting an error
    fx1 = -0.5*y1 + 0*x;
    fx2 = 4 - 0.3*y2 - 0.1*y1;
end