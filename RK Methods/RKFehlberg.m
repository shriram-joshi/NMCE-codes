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
   k2(i) = f(x(i) + h/5, y(i) + (k1(i)*h)/5);
   k3(i) = f(x(i) + (3*h)/4, y(i) + (3*k1(i)*h)/40 + (9*k2(i)*h)/40);
   k4(i) = f(x(i) + (3*h)/5, y(i) - (3*k1(i)*h)/10 - (9*k2(i)*h)/10 + (6*k3(i)*h)/5);
   k5(i) = f(x(i) + h, y(i) - (11*k1(i)*h)/54 + (5*k2(i)*h)/2 - (70*k3(i)*h)/27 + (35*k4(i)*h)/27);
   k6(i) = f(x(i) + (7*h)/8, y(i) + (1631*k1(i)*h)/55296 + (175*k2(i)*h)/512 + (575*k3(i)*h)/13824 + (44275*k4(i)*h)/110592 + (253/4096)*(k5(i)*h));
   % Fifth order estimate
   y(i+1) = y(i) + h*((2825/27648)*k1(i) + (18575/48384)*k3(i) + (13525/55296)*k4(i) + (277/14336)*k5(i) + (1/4)*k6(i));
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