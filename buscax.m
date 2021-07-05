function x=buscax(qi,a,b,m)
%buscax.m is used to evenly space nodes around mesh rows. Finds value of x
%that minimises the difference between the desired and true arc length at
%point x, by minimising function defined in ss1.m over x
%
%Inputs
%   qi = desired arc length up to point x
%   a,b = X- and Y- radii of mesh row 
%   m = exponent in equation (x/a)^m + y(/b)^m = 1 used to define shape of
%           mesh
%
%Outputs
%   x = x coordinate on ellipse that gives the required arc length
%

x = fzero(@ss1,[1e-6*a a-1e-6*a],[],0,qi,a,b,m);
                                %solution found as zero crossing of
                                %function defined in ss1.m

end