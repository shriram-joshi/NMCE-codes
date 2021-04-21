a=0;
b=0.8;

syms x
f = @(x) (0.2 + 25*x - 200*x.^2 + 675*x.^3 - 900*x.^4 + 400*x.^5);

avg = int(diff(diff(f,x)),[a b])/(b-a);
display(avg);