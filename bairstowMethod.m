clearvars;
%Inputs
pol = [4,12,13,12,9]; % in the form a0 + a1x + a2x^2 + ... + anx^n
es = 1;
r = -1.3333;
s = -0.4444;

%initialization 
iter = 0;
n = length(pol);
ea_r = 100;
ea_s = 100;

while (ea_r > es && ea_s > es)
    iter = iter + 1;
    disp("Iteration " + iter);
    
    % Finding b's and c's
    b(n) = pol(n);
    b(n-1) = pol(n-1) + r*b(n);
    c(n) = b(n);
    c(n-1) = b(n-1) + r*c(n);
    for i = n-2:-1:1
        b(i) = pol(i) + r*b(i+1) + s*b(i+2);
        c(i) = b(i) + r*c(i+1) + s*c(i+2);
    end
    
    disp("from b0 to bn");
    display(b); %from b0 to bn
    disp("from c0 to cn, note: c0 is not needed but still displayed");
    display(c); %from c0 to cn, note: c0 is not needed but still displayed
    
    det = c(3)*c(3) - c(4)*c(2);
    
    if(det ~= 0)
        dr = (-b(2)*c(3) + b(1)*c(4))/det;
        ds = (-b(1)*c(3) + b(2)*c(2))/det;
        disp("r(" + iter + ") =" + r + "+" + dr)
        disp("s(" + iter + ") =" + s + "+" + ds)
        r = r + dr;
        s = s + ds;
        if r~=0 
            ea_r = 100*abs(dr/r);
        end
        if s~=0 
            ea_s = 100*abs(ds/s);
        end
        disp("es_r(" + iter + ") =" + ea_r);
        disp("es_s(" + iter + ") =" + ea_s);
    else
        r = r+1;
        s = s+1;
        iter = 0;
    end 
end

disp("Final r =" + r);
disp("Final s =" + s);

coeff_x = [1,-r,-s];
x = roots(coeff_x);

display(x);


[q,rem] = deconv(flip(pol),coeff_x);
disp("q is in the order anx^n + ... + a1x + a0")
display(q);
display(rem);

