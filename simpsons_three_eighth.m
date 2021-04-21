clearvars;
clc;
disp("Simpson's 3/8th Rule");

%Enter basic info here
a=0;
b=0.8;
h = (b-a)/3;
x = zeros(1,n-1);

%Enter true value if mentioned else put 'trueValue = integral(@f,a,b);'
true_value = integral(@f,a,b);


for i = 1:1:4
    x(i)=a+(i-1).*h;
end

I = (3*h./8).*(f(x(1)) + 3.*f(x(2)) + 3.*f(x(3)) + f(x(4)));

Et = true_value-I;
perEt = 100*abs(Et/true_value);

%displaying results
disp("h = " + h)
disp("x = ");
disp(x);
disp("f(xi)=");
disp(f(x));
disp("True value = " + true_value);
disp("I = " + I);
disp("Et = " + Et);
disp("%Et = " + perEt);

function fx = f(x)
    %define your function here
    fx = (0.2 + 25*x - 200*x.^2 + 675*x.^3 - 900*x.^4 + 400*x.^5);
end


