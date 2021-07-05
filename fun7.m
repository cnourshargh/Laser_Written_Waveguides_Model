bcval=[];bcdof=[];
alpha=0.035;
betha=0.035
for i=1:length(gcoord)
    tita=atan(gcoord(i,2)/gcoord(i,1));
    ro=sqrt((a^2*b^2)/((cos(tita))^2*(b^2-a^2)+a^2)); %Radio de elipse en funcion de tita a y b
    r=(gcoord(i,2)^2+gcoord(i,1)^2)^0.5;
            if r<ro+0.00001*a 
            X=alpha*gcoord(i,1)
            Y=betha*gcoord(i,2)
            bcdof1=2*i-1
            bcdof2=2*i
            bcdof=[bcdof bcdof1 bcdof2] 
            bcval=[bcval X Y]
            else
            end
end