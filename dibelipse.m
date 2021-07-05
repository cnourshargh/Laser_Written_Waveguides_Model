p=dibelipse(a,b,n)
% programa que dibuja una elipse equispaciada el perimetro dividido n
% a= radio horizontal 
% b= radio vertical
% n= divisiones del perimetro.
p=zeros(n+1,2)
for i=1:n+1
    tita=(i-1)*pi/(4*n);
    r=sqrt(a^2*b^2/((cos(tita))^2*(b^2-a^2)+a^2));
    p(i,:)=[r*cos(tita) r*sin(tita)];
end

plot(p(:,1),p(:,2))



