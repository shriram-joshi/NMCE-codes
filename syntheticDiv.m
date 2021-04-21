%Inputs
pol = [-120,-46,79,-3,-7,1]; %Enter dividend in order a0 + a1x + a2x^2 ... + anx^n
n = length(pol);
t = 2; %Divisor of the form (x-t)

%Procedure
r = pol(n);
pol(n) = 0;
for i = n-1:-1:1
    s = pol(i);
    pol(i) = r;
    r = s + r*t;
end

display(r);
display(pol); %In the order b0 + b1x + a2x^2 ... + bn-1x^(n-1)