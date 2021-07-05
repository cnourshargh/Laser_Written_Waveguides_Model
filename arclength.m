function integrand = arclength(x,a,b,m)          
%arclength.m defines integrand for arc length calculation in fun1.m

    integrand= @(x,a,b,c) sqrt(1+(-x*b^2./(((1-x.^m/a^2)*b^m).^(1/m)*a^m)).^2);
    
end
