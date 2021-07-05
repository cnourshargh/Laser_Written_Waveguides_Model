function [nx,ny,nz,x,y,refractive_index,dx,Dnx,Dny,Dnz,photoelastic]=matfesa_adaptive_mesh(alpha,alpha2,betha,betha2,a,b,bad_element_boundary_y,ellipticity,positions,opt,rotation,track_index)
%----------------------------------------------------------------------------
%

%Flat ellipse deformations to simulate deformations and index change in crystal by
%fs laser interaction

%produces matlab error messages 6.5
warning off MATLAB:fzero:UndeterminedSyntax
warning off MATLAB:m_warning_end_without_block


q_ellipse_plot_size = 20;
mesh_res = 0.05;


[gcoord,nodes]=mesh1(a,b,q_ellipse_plot_size,mesh_res,bad_element_boundary_y,ellipticity);             %extracts node coordinates and connectivity from mesh
%meshplotter(gcoord,nodes);
[nel,xor]=size(nodes);                   % number of elements
nnel=9;                  % number of nodes per element
ndof=2;                  % number of dofs per node
nnode=length(gcoord);   % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element
emodule=0.2;             % elastic modulus [only useful for isotropic materials]
poisson=0.2;             % Poisson's ratio [only useful for isotropic materials]
nglx=3; ngly=3;          % 2x2 Gauss-Legendre quadrature
nglxy=nglx*ngly;         % number of sampling points per element
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
kinmtx2=zeros(3,edof);   % kinematic matrix
%matmtx=zeros(3,3);      % constitutive matrix
%----------------------------
%  force vector
%----------------------------
ff(1)=0;               % force applied at node 17 in y-axis
                      % force applied at node 18 in y-axis
%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------
[point2,weight2]=feglqd2(nglx,ngly);       % sampling points & weights
%matmtx=fematiso(4,emodule,poisson);        % compute constitutive matrix
[compliance, photoelastic, refractive_index] = material_selector(opt,rotation);

strainel=[];

for iel=1:nel           % loop for the total number of elements

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
x=point2(intx,1);                  % sampling point in x-axis
wtx=weight2(intx,1);               % weight in x-axis
for inty=1:ngly
y=point2(inty,2);                  % sampling point in y-axis
wty=weight2(inty,2);               % weight in y-axis
[shape,dhdr,dhds]=feisoq9(x,y);    % compute shape functions and
                                   % derivatives at sampling point
jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  % compute Jacobian
detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix
[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob); % derivatives w.r.t.
                                               % physical coordinate
kinmtx2=fekine2d(nnel,dhdx,dhdy);          % compute kinematic matrix
%------------------------------
%  compute element matrix
%------------------------------
k=k+kinmtx2'*compliance*kinmtx2*wtx*wty*detjacob;    % element matrix
end
end                                   % end of numerical integration loop
iel;                                   % K of each item is printed
k;     
index=feeldof(nd,nnel,ndof);% extract system dofs associated with element
kk=feasmbl1(kk,k,index);  % assemble element matrices 
end
%-----------------------------
%   apply boundary conditions
%  
%-----------------------------
[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);

%----------------------------
%  solve the matrix equation
%----------------------------
disp=kk\ff;

dx = 0.2; %mesh spacing for plot

[Txx,Tyy,Txy,x,y] = plotsurfrect(nodes,gcoord,disp,a,b,dx,q_ellipse_plot_size,ellipticity);

[Txx,Tyy,Txy,x,y] = plot_full_strain_map(Txx,Tyy,Txy,dx);

[Txx,Tyy,Txy] = position_ellipses(Txx,Tyy,Txy,dx,positions);

[Dnx,Dny,Dnz]=supind(Txx,Tyy,photoelastic,refractive_index,x,y,a,b,track_index,positions);

nx = Dnx+refractive_index(1,1);
ny = Dny+refractive_index(2,2);
nz = Dnz+refractive_index(3,3);







