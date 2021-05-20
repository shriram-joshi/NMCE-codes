clearvars;
clc;

format short;

x = [3,4.5,7,9];
fx = [2.5,1,2.5,0.5];

if (length(x) == length(fx))
    n = length(x);
    m = zeros(1,n);
    y = zeros(n,2);

    for i = 1:1:n-1
        m(i) = (fx(i+1) - fx(i))./(x(i+1)-x(i));
    end

    for i = 1:1:n
        y(i,1) = m(i);
        y(i,2) = fx(i)- m(i).*x(i);
    end

    disp("m=");
    disp(m);

    disp("y=");
    disp(y);
    plot(x,fx);
else
    disp("Your data is inconsistent. Length of fx and x array should be equal");
end