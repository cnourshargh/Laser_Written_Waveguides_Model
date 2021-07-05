function [Txx,Tyy,Txy,x,y]=plotsurfrect(nodes,gcoord,disp,a,b,dx,plot_size,m)
% This function plots strain maps 

xi=[0:dx:plot_size]; % X values of mapping 
yi=[0:dx:plot_size]; % Y values of mapping
[x,y]=meshgrid(xi,yi); %produces 2 meshes based on the values dictated above
[g,h]=size(x);
Txx=zeros(length(xi),length(yi)); %creates zeros matrix with dimensions of largest side of xi and yi
Tyy=Txx;
Txy=Txx;
z=0;
for i=1:length(xi)
    for j=1:length(yi)
        if (x(i,j)/a)^m+(y(i,j)/b)^m < 0.9999   %if coordinate s within damage elipse
            Txx(i,j)=0;                         %set strains to be zero
             Tyy(i,j)=0;
             Txy(i,j)=0;
        elseif (x(i,j)/a)^m+(y(i,j)/b)^m > 0.9999       %if coordinate is outside of ellipse
              z=z+1; 
              iel=fun8(nodes,gcoord,a,b,x(i,j),y(i,j));  %determines in which element the point xij yij lies
            %calculates coordinates of the corners of the ith element
            x1=gcoord(nodes(iel,1),1);
        	x2=gcoord(nodes(iel,2),1);
            x3=gcoord(nodes(iel,3),1);
            x4=gcoord(nodes(iel,4),1);
            y1=gcoord(nodes(iel,1),2);
            y2=gcoord(nodes(iel,2),2);
            y3=gcoord(nodes(iel,3),2);
            y4=gcoord(nodes(iel,4),2);    
           [ri,si]=fform(x1,x2,x3,x4,y1,y2,y3,y4,x(i,j),y(i,j)); %Evaluate isoparametric function form value at xij and yij
            ep=fun6(nodes,gcoord,disp,iel,ri,si); % Evaluate deformation at xij yij point
 %                 T(i,j)=iel;          
               Txx(i,j)=ep(1,1);   % Asign ex deformation
               Tyy(i,j)=ep(2,1);   % Asign ey deformation  
               Txy(i,j)=ep(3,1);   % Asign exy deformation
          end   
        
       end
end


%The following lines plot deformations

%figure;
%surf(x,y,Txx);
%figure;
%surf(x,y,Tyy);
%figure;
%surf(x,y,Txy);
end