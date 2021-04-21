clearvars;
clc;
format long e;

% The following function is for 3 points but can be expanded to use for any
% n. Just enter x and c in the following arrays
x = [0.7745966692414834, 0, -0.7745966692414834];
c = [0.5555555555555556, 0.8888888888888888, 0.5555555555555556];

%enter interval of intgration if any else set - a = -1 and b = 1
a = 0;
b = 10;

tildec = (b-a)/2*c;
tildex = (b-a)/2*x + (b+a)/2;

% define function here and write 'tildex' in place of x
f = (9.8*68.1/12.5)*(1-exp(-(12.5/68.1).*tildex));
%f = exp(tildex).*cos(tildex);
I = sum(tildec.*f);

disp("I = " + I);
