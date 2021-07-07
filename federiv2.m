function [dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob)
%federiv2.m determines the derivatives of the two dimensional
%isoparametric shape functions with respect to the physical coordinates.
%
%Inputs
%   nnel - number of nodes per element   
%   dhdr - derivative of shape functions w.r.t. natural coordinate r
%   dhds - derivative of shape functions w.r.t. natural coordinate s
%   invjacob - inverse of 2-D Jacobian matrix
%
%Outputs
%   dhdx - derivative of shape function w.r.t. physical coordinate x
%   dhdy - derivative of shape function w.r.t. physical coordinate y
%

%------------------------------
%  calculate shape function derivatives
%------------------------------

for i=1:nnel
 dhdx(i)=invjacob(1,1)*dhdr(i)+invjacob(1,2)*dhds(i);
 dhdy(i)=invjacob(2,1)*dhdr(i)+invjacob(2,2)*dhds(i);
 end
