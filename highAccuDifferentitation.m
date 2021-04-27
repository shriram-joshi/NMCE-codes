clearvars;
clc;
format long;

% Enter step size and x at which function is to be differentiatied
x_eval = 0.5;
h = 0.25;
x = [x_eval - 3*h, x_eval - 2*h, x_eval - h, x_eval, x_eval + h, x_eval + 2*h, x_eval + 3*h];

%Solving
i=4;

% Forward Finite Divided difference
f1F = (-f(x(i+2))+4*f(x(i+1))-3*f(x(i)))/(2*h);
f2F = (-f(x(i+3))+4*f(x(i+2))-5*f(x(i+1))+2*f(x(i)))/(h^2);
% Backward Finite Divided difference
f1B = (3*f(x(i))-4*f(x(i-1))+f(x(i-2)))/(2*h);
f2B = (2*f(x(i))-5*f(x(i-1))+4*f(x(i-2))-f(x(i-3)))/(h^2);
% Centered Finite Divided difference
f1C = (-f(x(i+2))+8*f(x(i+1))-8*f(x(i-1))+f(x(i-2)))/(12*h);
f2C = (-f(x(i+2))+16*f(x(i+1))-30*f(x(i))+16*f(x(i-1))-f(x(i-2)))/(12*(h^2));

% Displaying answer

display(x);

disp("Forward Finite Divided difference:");
disp("f'(" + x_eval +") = " + f1F);
disp("f''(" + x_eval +") = " + f2F);

disp("Backward Finite Divided difference:");
disp("f'(" + x_eval +") = " + f1B);
disp("f''(" + x_eval +") = " + f2B);

disp("Centered Finite Divided difference:");
disp("f'(" + x_eval +") = " + f1C);
disp("f''(" + x_eval +") = " + f2C);

function fx = f(x)
    %Enter your function here
    fx = (-0.1*x^4 - 0.15*x^3 - 0.5*x^2 -0.25*x + 1.2);
end

