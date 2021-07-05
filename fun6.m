function epsilon=fun6(nodes,gcoord,disp,iel,x,y)
%calculates deformations in in local x,y coordinates of element ie1
% iel=2; 
ndof=2;                         %number of degrees of freedom per node
nnel=9;                         %number of nodes per element

edof=nnel*ndof;                 %element degrees of freedom

for i=1:nnel
nd(i)=nodes(iel,i);             % extract connected node for (i*10)th element
xcoord(i)=gcoord(nd(i),1);      % extract x value of the node
ycoord(i)=gcoord(nd(i),2);      % extract y value of the node
end

[shape,dhdr,dhds]=feisoq9(x,y);    % compute shape functions and
                                   % derivatives at sampling point
                                   
jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  % compute Jacobian

detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix

[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob); % derivatives w.r.t.
                                               % physical coordinate

kinmtx2=fekine2d(nnel,dhdx,dhdy);           % kinematic matrix

index=feeldof(nd,nnel,ndof);% extract system dofs for the element

for i=1:edof
eldisp(i)=disp(index(i));
end

eps=kinmtx2*eldisp';
epsilon=eps;
% strains(:,w)=epsilon
end


