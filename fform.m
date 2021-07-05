function [ri,si]=fform(x1,x2,x3,x4,y1,y2,y3,y4,xi,yi)
% Function that calculates the function evaluated at point xi yi of a
% element in nodes x1 x2 x3 x4 y1 y2 y3 y4
A=[1 x1 y1 x1*y1;1 x2 y2 x2*y2;1 x3 y3 x3*y3;1 x4 y4 x4*y4];
Ainv=inv(A);
H1=[1 xi yi xi*yi]*Ainv(:,1);
H2=[1 xi yi xi*yi]*Ainv(:,2);
H3=[1 xi yi xi*yi]*Ainv(:,3);
H4=[1 xi yi xi*yi]*Ainv(:,4);
r1=-1;
r2=1;
r3=1;
r4=-1;
s1=-1;
s2=-1;
s3=1;
s4=1;
ri=H1*r1+H2*r2+H3*r3+H4*r4;
si=H1*s1+H2*s2+H3*s3+H4*s4;