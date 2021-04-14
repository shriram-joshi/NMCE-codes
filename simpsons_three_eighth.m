clearvars;
disp("Simpson's 3/8th Rulr");

%Enter basic info here
a=0;
b=0.8;
%Enter true value if mentioned else equate to = [];
true_value = 1.640533;

h = (b-a)/3;

for i = 1:1:4
    x(i)=a+(i-1).*h;
    %Enter your function here
    fx(i) = (0.2 + 25*x(i) - 200*x(i).^2 + 675*x(i).^3 - 900*x(i).^4 + 400*x(i).^5);
end

I = (3*h./8).*(fx(1) + 3.*fx(2) + 3.*fx(3) + fx(4));

if ~isempty(true_value)
    Et = true_value-I;
    perEt = 100*abs(Et/true_value);
    disp("Et = " + Et);
    disp("%Et = " + perEt);
end

display(x);
display(fx);
disp("I = " + I);


