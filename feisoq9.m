function [shapeq9,dhdrq9,dhdsq9]=feisoq9(r,s)
%feisoq9.m computes the isoparametric 9 node quadrilateral shape functions
%at the selected (integration) point in terms of the natural coordinates.
%
%Inputs
%   r - X coordinate value of the selected point   
%   s - Y coordinate value of the selected point

%  Variable Description:
%     shapeq9 - shape functions for nine-node element
%     dhdrq9 - derivatives of the shape functions w.r.t. r
%     dhdsq9 - derivatives of the shape functions w.r.t. s

%
%  Notes:
%     1st node at (-1,-1), 2nd node at (1,-1)
%     3rd node at (1,1), 4th node at (-1,1)
%------------------------------------------------------------------------

% shape functions

 shapeq9(1)=0.25*(r^2-r)*(s^2-s);
 shapeq9(2)=0.25*(r^2+r)*(s^2-s);
 shapeq9(3)=0.25*(r^2+r)*(s^2+s);
 shapeq9(4)=0.25*(r^2-r)*(s^2+s);
 shapeq9(5)=0.5*(1-r^2)*(s^2-s);
 shapeq9(6)=0.5*(r^2+r)*(1-s^2);
 shapeq9(7)=0.5*(1-r^2)*(s^2+s);
 shapeq9(8)=0.5*(r^2-r)*(1-s^2);
 shapeq9(9)=(1-r^2)*(1-s^2);
 
 % derivatives

 dhdrq9(1)=0.25*(2*r-1)*(s^2-s);
 dhdrq9(2)=0.25*(2*r+1)*(s^2-s);
 dhdrq9(3)=0.25*(2*r+1)*(s^2+s);
 dhdrq9(4)=0.25*(2*r-1)*(s^2+s);
 dhdrq9(5)=0.5*(-2*r)*(s^2-s);
 dhdrq9(6)=0.5*(2*r+1)*(1-s^2);
 dhdrq9(7)=0.5*(-2*r)*(s^2+s);
 dhdrq9(8)=0.5*(2*r-1)*(1-s^2);
 dhdrq9(9)=(-2*r)*(1-s^2);
 
 
 dhdsq9(1)=0.25*(r^2-r)*(2*s-1);
 dhdsq9(2)=0.25*(r^2+r)*(2*s-1);
 dhdsq9(3)=0.25*(r^2+r)*(2*s+1);
 dhdsq9(4)=0.25*(r^2-r)*(2*s+1);
 dhdsq9(5)=0.5*(1-r^2)*(2*s-1);
 dhdsq9(6)=0.5*(r^2+r)*(-2*s);
 dhdsq9(7)=0.5*(1-r^2)*(2*s+1);
 dhdsq9(8)=0.5*(r^2-r)*(-2*s);
 dhdsq9(9)=(1-r^2)*(-2*s);
 
 
