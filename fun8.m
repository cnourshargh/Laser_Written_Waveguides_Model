function elemento=fun8(nodes,gcoord,a,b,xij,yij)
%function to determine in which element the point xij yij lies
[m,n]=size(nodes);
for ele=1:m
    hi(ele)=mean([((gcoord(nodes(ele,9),1)-xij).^2+(gcoord(nodes(ele,9),2)-yij).^2).^0.5;
                  ((gcoord(nodes(ele,1),1)-xij).^2+(gcoord(nodes(ele,1),2)-yij).^2).^0.5; 
                  ((gcoord(nodes(ele,2),1)-xij).^2+(gcoord(nodes(ele,2),2)-yij).^2).^0.5;
                  ((gcoord(nodes(ele,3),1)-xij).^2+(gcoord(nodes(ele,3),2)-yij).^2).^0.5;
                  ((gcoord(nodes(ele,4),1)-xij).^2+(gcoord(nodes(ele,4),2)-yij).^2).^0.5;
                  ((gcoord(nodes(ele,5),1)-xij).^2+(gcoord(nodes(ele,5),2)-yij).^2).^0.5
                  ((gcoord(nodes(ele,6),1)-xij).^2+(gcoord(nodes(ele,6),2)-yij).^2).^0.5;
                  ((gcoord(nodes(ele,7),1)-xij).^2+(gcoord(nodes(ele,7),2)-yij).^2).^0.5;
                  ((gcoord(nodes(ele,8),1)-xij).^2+(gcoord(nodes(ele,8),2)-yij).^2).^0.5]);

end
[j,elemento]=min(hi);