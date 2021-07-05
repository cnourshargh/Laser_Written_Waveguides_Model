%fun 4
%programa para mostrar strainxx en funcion de x para y=0
def=zeros(n2+1,3)
def(1,:)=strainel(1,:)
def(n2+1,:)=strainel(length(strainel),:)
e=2;
for i=3:2:length(strainel)-1
        def(e,:)=strainel(i,:)+strainel(i-1,:)
        e=e+1;
end
figure
plot(def(:,1))
figure
plot(def(:,2),'r')  %Ploteo strain y
figure
plot(def(:,3),'g')   %Ploteo strain xy



