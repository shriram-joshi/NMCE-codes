format long e;
clearvars;
toll=0.000001;
a=0;
b=0.8;
c=(b+a)/2;
I = qStep(a,c,b,toll);

%displaying results
disp("Real value = " + integral(@f,a,b));
display(I);

function eval = qStep(a,c,b,toll)
    h1 = b-a;
    h2 = h1/2;
    d = (a+c)/2;
    e = (c+b)/2;
    I1 = (h1/6)*(f(a)+4*f(c)+f(b));
    I2 = (h2/6)*(f(a)+4*f(d)+2*f(c)+4*f(e)+f(b));
    
    if(abs(I2-I1) <= toll)
        I = I2 + (I2-I1)/15;
    else 
        Ia = qStep(a,d,c,toll);
        Ib = qStep(c,e,b,toll);
        I = Ia+Ib;        
    end
    eval = I;
end

function fx = f(x)
    fx = (0.2 + 25*x - 200*x.^2 + 675*x.^3 - 900*x.^4 + 400*x.^5);
end

