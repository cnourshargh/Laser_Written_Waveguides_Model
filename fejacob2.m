function [jacob2]=fejacob2(nnel,dhdr,dhds,xcoord,ycoord)
%fejacob2.m produces the Jacobian for the two dimensional mapping between
%the finite element coordinate system and the cartesian coordiante system
%
%Inputs
%   nnel - number of nodes per element   
%   dhdr - derivative of shape functions w.r.t. natural coordinate r
%   dhds - derivative of shape functions w.r.t. natural coordinate s
%   xcoord - x axis coordinate values of nodes
%   ycoord - y axis coordinate values of nodes
%
%Outputs
%   jacob2 - Jacobian for one-dimension
%   

%--------------------------------
%  Initialise Jacobian matrix
%--------------------------------

jacob2=zeros(2,2);                         

%--------------------------------
%  define terms in matrix
%--------------------------------

 for i=1:nnel
 jacob2(1,1)=jacob2(1,1)+dhdr(i)*xcoord(i);
 jacob2(1,2)=jacob2(1,2)+dhdr(i)*ycoord(i);
 jacob2(2,1)=jacob2(2,1)+dhds(i)*xcoord(i);
 jacob2(2,2)=jacob2(2,2)+dhds(i)*ycoord(i);
 end
