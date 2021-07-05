%DIbuja la deformada.
m=0;
for i=1:2:length(disp)
    m=m+1;
    dispx(m)=disp(i);
end
m=0;
for i=2:2:length(disp)
    m=m+1;
    dispy(m)=disp(i);
end
figure
plot(gcoord(:,1),gcoord(:,2),'o',gcoord(:,1)+dispx',gcoord(:,2)+dispy','+');


