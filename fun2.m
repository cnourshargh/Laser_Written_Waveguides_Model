function [bcdof,bcval]=fun2(gcoord,nodes,a,b,alpha,alpha2,betha,betha2,m)
%fun2.m defines the boundary strains and constraints of all the
%nodes in the finite element mesh to force symmetry constraints. This
%allows consideration of only one quadrant (a quarter ellipse). Each node
%is indexed by two numbers, one relating to its movement in x, and one to
%its movement in y.
%
%Inputs
%   gcoord = coordinates of all nodes in finite element mesh
%   nodes = lists nodes contained in each finite element
%   a,b = X- and Y- radii of ellipse boundary
%   alpha, betha = linear expansion coefficient applied to strain (track)
%                   boundary; can be set differently for X- and Y- but not
%                   done here
%   alpha2,betha2 = quadratic expansion coefficient applied to strain (track)
%                   boundary; can be set differently for X- and Y- but not
%                   done here
%   m = ellipticity; n in equation (x/a)^n + y(/b)^n = 1 used to define strain
%           boundary
%
%Outputs
%   bcdof = vector containing index of all node degrees of freedom which
%           are constrained
%   bcval = vector containing strain of each of the nodes identified in
%           bcdof
%

w=1;                            %counter used to position values in bc vectors

for i=1:length(gcoord)          %for every node in the finite element mesh
    if gcoord(i,2)<0.0001 && (((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m < 0.99999) || (((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m) > 1.00001))
                                %if node lies on X axis and is not on
                                %elliptical strain boundary
                                
            bcdof(w)=2*i;       %index for Y movement of node i
            bcval(w)=0;         %zero strain in y, hence still lies on X axis
            w=w+1;              %increments counter by 1
            
    elseif gcoord(i,1)<0.0001 && (((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m < 0.99999) || (((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m) > 1.00001))
                                %if node lies on Y axis and is not on
                                %elliptical strain boundary
                                
            bcdof(w)=2*i-1;     %index for X movement of node i      
            bcval(w)=0;         %zero strain in X, hence still lies on Y axis
            w=w+1;              %increments counter by 1
            
    elseif ((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m > 0.9999) && (((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m) < 1.0001)
                                %if node lies on elliptical strain boundary
                                
            bcdof(w)=2*i-1;     %index for X movement of node i
			bcval(w)=alpha*gcoord(i,1)+alpha2*gcoord(i,1)^2;
                                %strain in x dictated by horizontal expansion
                                %parameters
                                
            bcdof(w+1)=2*i;     %index for movement in y
            bcval(w+1)=betha*gcoord(i,2)+betha2*gcoord(i,2)^2;
                                %strain in Y dictated by vertical expansion
                                %parameters
                                
            w=w+2;              %counter increments by 2
    end 
end





