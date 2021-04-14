clearvars;
clc;
disp("Simpson's 1/3rd Rule");

%Enter basic info here
a=0;
b=0.8;
n=4;
h = (b-a)/n;
disp("h = " + h)
%Enter true value if mentioned else equate to = [];
true_value = 1.640533;

for i = 1:1:n+1
    x(i)=a+(i-1).*h;
    %Enter your function here
    fx(i) = (0.2 + 25*x(i) - 200*x(i).^2 + 675*x(i).^3 - 900*x(i).^4 + 400*x(i).^5);
end

I = (h./3).*(fx(1) + fx(n+1));

for i = 2:1:n
    if rem((i-1),2)==0
        I = I + (2*h/3)*fx(i);
    else 
        I = I + (4*h/3)*fx(i);
    end
end

if ~isempty(true_value)
    Et = true_value-I;
    perEt = 100*abs(Et/true_value);
    disp("Et = " + Et);
    disp("%Et = " + perEt);
end

display(x);
display(fx);
disp("I = " + I);


