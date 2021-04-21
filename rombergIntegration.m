clearvars;
clc;
format long;

a=0;
b=0.8;
n=1;
maxit = 3;
Es = 0.00001;
Ea = 100;
I = zeros(n+1,n+1);
iter = 0;

I(1,1) = trapInt(n,a,b);

while iter <= maxit-1 && Ea >= Es 
    iter = iter +1;
    n = 2^iter;
    I(iter+1,1) = trapInt(n,a,b);
    for k = 2:iter+1
        j = 2+iter-k;
        I(j,k) = ((4^(k-1))*I(j+1,k-1)-I(j,k-1))/((4^(k-1))-1);
    end
    Ea = abs((I(1,iter+1)-I(2,iter))/I(1,iter+1))*100;
end

disp("I =");
disp(I);

I_romberg = I(1,iter+1);

disp("I(O(h^"+ (4+2*(maxit-1)) + ")) = " + I_romberg);
display(Ea);


function tIntegral = trapInt(n,a,b)
    h = (b-a)/n;
    x = zeros(1,n-1);

    for i = 1:1:n+1
        x(i)=a+(i-1).*h;
    end

    tIntegral = (h./2).*(f(x(1)) + f(x(n+1)));

    for i = 2:1:n
        tIntegral = tIntegral + h.*f(x(i));
    end
end

function fx = f(x)
    %define your function here
    fx = double(0.2 + 25*x - 200*x^2 + 675*x^3 - 900*x^4 + 400*x^5);
end