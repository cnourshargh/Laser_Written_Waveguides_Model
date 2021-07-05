function gc=fun1(a,b,bad_element_boundaries,b0,m)
%fun1.m calculates the coordinates of the nodes on the eleiptical mesh. It
%is used for both the simple and the adaptive meshing so has the capability
%for both evenly spacing the nodes and for adding resolution in a narrow
%region.
%
%Inputs
%   a,b = X- and Y- radii of mesh row
%   bad_element_boundaries = used in matfesa_adaptive_mesh.m. Coordinates 
%                               on the strain boundary between which a
%                               higher resolution is required. left empty
%                               for matfesa_simple_mesh.m.
%   b0 = Y- radius of boundary ellipse, used for calcualting region in
%           which adaptive mesh is needed
%   m = ellipticity in equation (x/a)^m + y(/b)^m = 1 used to define shape of finite 
%           mesh row     
%                               
%Outputs
%   gc = collumn vector of node coordinates
%

%-----------------------------------------
% Calculate mesh regions and define node densities
%-----------------------------------------

if isempty(bad_element_boundaries)          %if bad element boundary isn't specified
                                            %ie simple mesh
                                            
    diva=8;                                 %number of values in which section a (lower) is divided 
    divb=8;                                 %number of values in which section b (middle) is divided
    divc=8;                                 %number of values in which section c (upper) is divided

    h1=0.5*b;                               %height at which section a ends and section b begins
    h2=0.9*b;                               %height at which section b ends and section c begins

    
elseif bad_element_boundaries(2)==0         %if bad element boundary is against the X- axis
    diva = 40;                              %number of values in which section a (lower) is divided 
    divb = 8;                               %number of values in which section b (middle) is divided 
    divc = 8;                               %number of values in which section c (upper) is divided 
    
    h1 = (bad_element_boundaries(4)/b0)*b;  %section a ends at the top of the bad element 
    h2 = 0.9*b;                             %section c starts starts as in simple mesh
        
elseif bad_element_boundaries(3)==0         %if bad element boundary is against Y- axis
    diva = 8;                               %number of values in which section a (lower) is divided 
    divb = 8;                               %number of values in which section b (middle) is divided 
    divc = 40;                              %number of values in which section c (upper) is divided 
    
    h1 = 0.5*b;                             %section a ends as in simple mesh
    h2 = (bad_element_boundaries(2)/b0)*b;  %section c begins at the beginning of the bad element 

else                                        %if bad element boundary is not against either axis
    
    diva=8;                                 %number of values in which section a (lower) is divided 
    divb=40;                                %number of values in which section b (middle) is divided 
    divc=8;                                 %number of values in which section c (upper) is divided 

    h1=(bad_element_boundaries(2)/b0)*b;	%section b starts at the bottom of the bad element
    h2=(bad_element_boundaries(4)/b0)*b; 	%section c starts at the top of the bad element
end

%-----------------------------------------
% Calculate node coordinates
%-----------------------------------------

x1=((1-(h1^m)/b^m)*a^m)^(1/m);              %calculates x coordinates of point given by h1 
x2=((1-(h2^m)/b^m)*a^m)^(1/m);              %calculates x coordinates of point given by h2

integrand= @(x,a,b,m) sqrt(1+(-x*b^m./(((1-x.^m/a^2)*b^m).^(1/m)*a^m)).^2);
                                            %equation to be integrated to
                                            %find arc length

QA = integral(@(x) integrand(x,a,b,m),x1,a-a*10^-9);   
                                            %solves integral over section a
qi=QA/diva;                                 %length of each sub arc

QB = integral(@(x) integrand(x,a,b,m),x2,x1);  
                                            %solves integral over section b
qj=QB/divb;                                 %length of each sub arc


QC = integral(@(x) integrand(x,a,b,m),0,x2);  
                                            %solves integral over section c
qk=QC/divc;                                 %length of each sub arc

Qtot=QA+QB+QC;                              %calculates total arc length

xr(1)=a;                                    %sets coordinate of radius on X- axis
yr(1)=0;
t=2;                                        %counter to be used in following itteration

for g=1:diva                                %itterates over nodes in section a
    qr=Qtot-g*qi;                           %calculates arclength upto node g
    xr(t)=buscax(qr,a,b,m);                 %calulates X- coordinate of node
    yr(t)=((1-xr(t)^m/a^m)*b^m)^(1/m);      %calulates Y- coordinate of node
    t=t+1;                                  %itterates counter
end

for g=1:divb                                %itterates over nodes in section b
    qr=Qtot-(diva*qi+g*qj);                 %calculates arclength upto node g
    xr(t)=buscax(qr,a,b,m);                 %calulates X- coordinate of node
    yr(t)=((1-xr(t)^m/a^m)*b^m)^(1/m);      %calulates Y- coordinate of node
    t=t+1;                                  %itterates counter
end

for g=1:(divc-1)                            %itterates over nodes in section b
    qr=Qtot-(diva*qi+divb*qj+g*qk);         %calculates arclength upto node g
    xr(t)=buscax(qr,a,b,m);                 %calulates X- coordinate of node
    yr(t)=((1-xr(t)^m/a^m)*b^m)^(1/m);      %calulates Y- coordinate of node
    t=t+1;                                  %itterates counter
end

xr(t)=0;                                    %sets coordinate of radius on Y- axis
yr(t)=b;

gc=[xr' yr'];                               %saves node coordinates as collumn vectors in gc



