function [kinmtx2]=fekine2d(nnel,dhdx,dhdy)
% fekine2d.m determines the coefficients of the kinemtaic equation relating the strans and
% displacemts for the current element.
%
%Inputs
%   nnel - number of nodes per element
%   dhdx - derivatives of shape functions with respect to x   
%   dhdy - derivatives of shape functions with respect to y
%
%Outputs
%   kinmtx2 - coefficients for the kinematic equation for the current
%               element
%
 for i=1:nnel
 i1=(i-1)*2+1;  
 i2=i1+1;
 kinmtx2(1,i1)=dhdx(i);
 kinmtx2(2,i2)=dhdy(i);
 kinmtx2(3,i1)=dhdy(i);
 kinmtx2(3,i2)=dhdx(i);
 end
