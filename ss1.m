function s=ss1(xf,xi,qi,a,b,m) 
%ss1.m is used to evenly space the nodes in the finite element mesh along
%each row. used by buscax.m in fun1.m
%
%Inputs
%   xf,xi = final and initial values of x, between which integral is
%           computed
%   qi = arc length at point being calculated
%   a,b = X- and Y- radii of mesh row 
%   m = exponent in equation (x/a)^m + y(/b)^m = 1 used to define shape of
%           mesh
% 
%Outputs
%   s = difference between desired arc length and arc length at position x

integrand= @(x,a,b,c) sqrt(1+(-x*b^m./(((1-x.^m/a^2)*b^m).^(1/m)*a^m)).^2);
                            %function to be integrated to find arc length
                            %at point x
s = qi- integral(@(x) integrand(x,a,b,m),xi,xf);
                            %difference between desired arc length and arc
                            %length at point x
end
