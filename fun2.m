function [bcdof,bcval]=fun2(gcoord,nodes,a,b,alpha,alpha2,betha,betha2,m)
%program calculates symmetry constraints
w=1;
for i=1:length(gcoord)
    if gcoord(i,2)<0.0001 && (((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m < 0.99999) || (((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m) > 1.00001))
            bcdof(w)=2*i;
            bcval(w)=0;
            w=w+1;
    elseif gcoord(i,1)<0.0001 && (((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m < 0.99999) || (((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m) > 1.00001))
            bcdof(w)=2*i-1;             %the same for horizontal displacements at y = 0
            bcval(w)=0;
            w=w+1;
    elseif ((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m > 0.9999) && (((gcoord(i,1)/a)^m+(gcoord(i,2)/b)^m) < 1.0001)
            bcdof(w)=2*i-1;
			bcval(w)=alpha*gcoord(i,1)+alpha2*gcoord(i,1)^2;
            bcdof(w+1)=2*i;
            bcval(w+1)=betha*gcoord(i,2)+betha2*gcoord(i,2)^2;
            w=w+2;     
    end 
end





