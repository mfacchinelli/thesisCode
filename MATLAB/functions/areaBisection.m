function c = areaBisection(t,acc,a0,b0)

%...Function handles
area_find = @(a) [trapz(t(1:a),acc(1:a)),trapz(t(a+1:end),acc(a+1:end))];
area_comp = @(a) max(abs(a(1)/a(2)),abs(a(2)/a(1)));

%...Iterative scheme
iter = 0;
while iter <= 100
    iter = iter+1;
    a = a0;
    b = b0;
    
    area1 = area_find(a);
    area2 = area_find(b);
    
    ratio1 = area_comp(area1);
    ratio2 = area_comp(area2);
    
    c = round((a + b)/2);
    if ratio1 > ratio2
        d = b;
        a0 = c;
    else
        d = a;
        b0 = c;
    end
    
    if abs(d-c) <= 1
        fprintf('Converged after %d iterations.\n',iter)
        break
    end
end