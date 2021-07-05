function [matmtrx]=fematiso(iopt,elastic,poisson)

%------------------------------------------------------------------------
%  Purpose:
%     determine the constitutive equation for isotropic material
%
%  Synopsis:
%     [matmtrx]=fematiso(iopt,elastic,poisson) 
%
%  Variable Description:
%     elastic - elastic modulus
%     poisson - Poisson's ratio   
%     iopt=1 - plane stress analysis
%     iopt=2 - plane strain analysis
%     iopt=3 - axisymmetric analysis
%     iopt=4 - three dimensional analysis 
%------------------------------------------------------------------------

 if iopt==1        % plane stress
   matmtrx= elastic/(1-poisson*poisson)* ...
   [1  poisson 0; ...
   poisson  1  0; ...
   0  0  (1-poisson)/2];

 elseif iopt==2     % plane strain
   matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson 0; 
   poisson  (1-poisson)  0;
   0  0  (1-2*poisson)/2];

 elseif iopt==3     % axisymmetry
   matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson  0; 
   poisson  (1-poisson)   poisson  0;
   poisson  poisson  (1-poisson)   0;
   0    0    0   (1-2*poisson)/2];
 
 elseif iopt==4 %Lithium Niobate orthotropic approx. Cut x, x horizontal, z-axis vertical and y-axis normal to plane
 % Approx. Orthotropic of Lithium Niobate z-cut, horizontal x, z-axis
 % vertical and y-axis normal to plane
 %matmtrx=1e-1*[ 2.0449 0.685 0;0.685 2.4675 0;0 0 0.609];
  
 %Lithium Niobate orthotropic Approx x cut, horizontal z, z axis
 % horizontal and y-axis and normal to plane 
 matmtrx= [0.2446    0.071   0;0.071   0.2058     0;0     0    0.0602];

 
 else     % three-dimension
   matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson   0   0    0; 
   poisson  (1-poisson)   poisson   0   0    0;
   poisson  poisson  (1-poisson)    0   0    0;
   0    0    0    (1-2*poisson)/2   0    0;
   0    0    0    0    (1-2*poisson)/2   0;
   0    0    0    0    0   (1-2*poisson)/2];

 end
