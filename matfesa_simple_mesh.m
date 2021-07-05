function [bad_element_boundary]=matfesa_simple_mesh(alpha,alpha2,betha,betha2,a,b,ellipticity,opt,rotation)
%matfesa_simple_mesh.m is used as the first stage in the adaptive meshing
%procedure. It identifies the the location with the greatest error produced
%by the finite element photoelasticity code. matfesa_adaptive_mesh.m is
%then run with many more elements where ther error ocurs to minimise
%erroneous area so the error can be removed with a median filter
%
%Inputs
%   alpha, betha = linear expansion coefficient applied to strain (track)
%                   boundary; can be set differently for X- and Y- but not
%                   done here
%   alpha2, betha2 = quadratic expansion coefficient applied to strain (track)
%                   boundary; can be set differently for X- and Y- but not
%                   done here
%   a,b = track X- and Y- radii respectively
%   ellipticity = n in equation (x/a)^n + y(/b)^n = 1 used to define strain
%               boundary
%   opt = Used to select material parameters, see material_selector.m for
%       more
%   rotation = used to define material orientation when not specified in
%               material_selector.m
%Outputs
%   bad_element_boundary = X- and Y- coordinates of pioints on elliptical
%                           boundary between which finite elements must be
%                           added by adaptive mesh

%produces matlab error messages 6.5
warning off MATLAB:fzero:UndeterminedSyntax
warning off MATLAB:m_warning_end_without_block


%-----------------------------------------
%define mesh dimesions
%-----------------------------------------
q_ellipse_plot_size = b+3;  %used to determine area over which elements are distributed
mesh_res = 0.15;            %used to determine ellipse row spacing

%-----------------------------------------
%generate mesh
%-----------------------------------------

[gcoord,nodes]=mesh1(a,b,q_ellipse_plot_size,mesh_res,[],ellipticity);         
                                        %extracts node coordinates and connectivity from mesh 
meshplotter(gcoord,nodes);             %plots mesh in X and Y

%-----------------------------------------
%define program parameters
%-----------------------------------------
[nel,xor]=size(nodes);                  % number of elements
nnel=9;                                 % number of nodes per element
ndof=2;                                 % number of dofs per node
nnode=length(gcoord);                   % total number of nodes in system
sdof=nnode*ndof;                        % total system dofs  
edof=nnel*ndof;                         % degrees of freedom per element
emodule=0.2;                            % elastic modulus [only useful for isotropic materials]
poisson=0.2;                            % Poisson's ratio [only useful for isotropic materials]
nglx=3; ngly=3;                         % 2x2 Gauss-Legendre quadrature
nglxy=nglx*ngly;                        % number of sampling points per element

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------
[bcdof,bcval]=fun2(gcoord,nodes,a,b,alpha,alpha2,betha,betha2,ellipticity);

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------
ff=zeros(sdof,1);       % system force vector
kk=zeros(sdof,sdof);    % system matrix
disp=zeros(sdof,1);     % system displacement vector
eldisp=zeros(edof,1);   % element displacement vector
stress=zeros(nglxy,3);  % matrix containing stress components
strain=zeros(nglxy,3);  % matrix containing strain components
index=zeros(edof,1);    % index vector
kinmtx2=zeros(3,edof);  % kinematic matrix

%----------------------------
%  force vector
%----------------------------
ff(1)=0;               % force applied at node 17 in y-axis
                       % force applied at node 18 in y-axis
                      
%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------
[point2,weight2]=feglqd2(nglx,ngly);       % sampling points & weights

[compliance, photoelastic, refractive_index] = material_selector(opt,rotation);
                            %uses material selector code to form 2D
                            %material matrices
strainel=[];

for iel=1:nel               % loop for the total number of elements

for i=1:nnel
nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
end

k=zeros(edof,edof);         % initialization of element matrix to zero

%--------------------------------
%  numerical integration
%--------------------------------
for intx=1:nglx
x=point2(intx,1);                               % sampling point in x-axis
wtx=weight2(intx,1);                            % weight in x-axis
for inty=1:ngly
y=point2(inty,2);                               % sampling point in y-axis
wty=weight2(inty,2);                            % weight in y-axis
[shape,dhdr,dhds]=feisoq9(x,y);                 % compute shape functions and
                                                % derivatives at sampling point
jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  % compute Jacobian
detjacob=det(jacob2);                           % determinant of Jacobian
invjacob=inv(jacob2);                           % inverse of Jacobian matrix
[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob);  % derivatives w.r.t.
                                                % physical coordinate
kinmtx2=fekine2d(nnel,dhdx,dhdy);               % compute kinematic matrix

%------------------------------
%  compute element matrix
%------------------------------
k=k+kinmtx2'*compliance*kinmtx2*wtx*wty*detjacob;   % element matrix
end
end                                             % end of numerical integration loop
   
index=feeldof(nd,nnel,ndof);                    % extract system dofs associated with element
kk=feasmbl1(kk,k,index);                        % assemble element matrices 
end

%-----------------------------
%   apply boundary conditions
%-----------------------------
[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);            %applies strain boundary conditions by zeroing some elements

%----------------------------
%  solve the matrix equation
%----------------------------
disp=kk\ff;                                     %solves for displacements of nodes                             

%--------------------------
%solve elasto-optic problem
%--------------------------
dx = 0.1;                                       %mesh spacing for plot

[Txx,Tyy,Txy,x,y]=plotsurfrect(nodes,gcoord,disp,a,b,dx,q_ellipse_plot_size,ellipticity);
                                                %Calculates strain for each point in cartesian mesh 

[Dnx,Dny,Dnz]=supind(Txx,Tyy,photoelastic,refractive_index,x,y,a,b);              
                                                %solves for change in refractive index

%-----------------------------------------
%find coordinates of error region
%-----------------------------------------
n_elements_per_row = 12;                        %defines number of elements in each row                 

A = real(Dnx)+real(Dny);                        %sums refractive index fields

ddx_A = (1/dx)*(A(1:end-1,:)-A(2:end,:));       %finds derivative of sum in x to find slope

ddy_A = (1/dx)*(A(:,1:end-1)-A(:,2:end));       %finds derivative of sum in x to find slope

index = find(abs(ddx_A) == max(ddx_A(:)));      %finds location of max in x derivative

worst_element = fun8(nodes,gcoord,a,b,x(index),y(index));
                                                %finds element containing max

badelement = mod(worst_element,n_elements_per_row);
                                                %finds corresponding element on shape circumference

if badelement == 0                              %if bad element is a multiple of the number of elements per row
                                                %so is against y axis
    badelement = n_elements_per_row;            %elements index from 1 so to n so moves 0 to n
    bad_element_boundary = [gcoord(nodes(badelement,1)-4,1),gcoord(nodes(badelement,1)-4,2),gcoord(nodes(badelement,4),1),gcoord(nodes(badelement,4),2)];
                                                %sets boundary from bottom of element n-2 to top of element n
                                                
elseif badelement == 1                          %if bad element is against x axis         
    bad_element_boundary = [gcoord(nodes(badelement,1),1),gcoord(nodes(badelement,1),2),gcoord(nodes(badelement,4)+4,1),gcoord(nodes(badelement,4)+4,2)];
                                                %sets boundary from bottom of element n to top of element n+2                                            

else
    bad_element_boundary = [gcoord(nodes(badelement,1)-2,1),gcoord(nodes(badelement,1)-2,2),gcoord(nodes(badelement,4)+2,1),gcoord(nodes(badelement,4)+2,2)];
                                                %sets boundary from bottom of element n-1 to top of element n+1
end






