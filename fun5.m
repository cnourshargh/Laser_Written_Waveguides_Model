function epsilon=fun5(edof,ndof,nnel,iel,x,y)

%Pos procesado para calcular las deformaciones en la coordenada local x,y
%del elemento iel. 
%iel=8; 
%edof=8;
% ndof=2;
% nnel=4;
% x=0.5774;  % coordenada local de punto donde quiero averiguar las tensiones o deformaciones
% y=0.5774;

for i=1:nnel
nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
end

[shape,dhdr,dhds]=feisoq4(x,y);    % compute shape functions and
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

epsilon=kinmtx2*eldisp;

