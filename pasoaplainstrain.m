%reduced matrix calculation to show plane deformations
s=inv(Cab);
alpha1=(s(4,1)-(s(1,3)*s(3,1))/s(3,3))/(-s(4,4)+(s(3,2)*s(1,3))/s(3,3));
alpha2=(s(4,2)-(s(1,3)*s(3,1))/s(3,3))/(-s(4,4)+(s(3,2)*s(1,3))/s(3,3));
Sp11=(s(1,1)-(s(1,3)*s(1,3))/s(3,3))+((s(1,4)-s(1,3)*s(3,4))/s(3,3))*alpha1;
Sp12=(s(1,2)-(s(3,2)*s(1,3))/s(3,3))+((s(1,4)-s(1,3)*s(3,4))/s(3,3))*alpha2;
Sp21=(s(2,1)-(s(1,3)*s(1,3))/s(3,3))+((s(1,4)-s(1,3)*s(3,4))/s(3,3))*alpha1;
Sp22=(s(2,2)-(s(3,2)*s(1,3))/s(3,3))+((s(2,4)-s(1,3)*s(3,4))/s(3,3))*alpha2;
Sp33=s(6,6)-(s(6,5)*s(5,6))/s(5,5);

Sp=[ Sp11 Sp12 0; Sp21 Sp22 0; 0 0 Sp33]
Cp=inv(Sp)