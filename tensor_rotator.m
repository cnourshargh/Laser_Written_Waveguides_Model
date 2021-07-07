function new_mtx = tensor_rotator(B,rotation)
%tensor_rotator.m is used to rotate the stiffness and photelastic matrices in 
%material_selector.m. It constructs 3 rotation matrices, relating to the
%three axes about which rotations are defined and multiplies all three
%together. Values from this rotation matrix are extracted to produce a
%matrix that can be used to rotate the 6x6 matrices.
%
%Inputs
%   B = 6x6 matrix to be rotated
%   rotation = vector describing matrix rotation, of form [theta, phi, alpha].       
%       theta = rotation of material matrices about the x axis
%       phi = rotation of material matrices about the y axis
%       alpha = rotation of material matrices about the z axis
%
%Outputs
%   new_mtx = B rotated into new frame
%

%--------------------------------
%  identify rotation angles from input
%--------------------------------

theta = rotation(1);
phi = rotation(2);
alpha = rotation(3);

%--------------------------------
%  construct separate rotation matrices
%--------------------------------

R_x = [1, 0, 0;
    0, cosd(theta), -sind(theta);
    0, -sind(theta), cosd(theta)];

R_y = [cosd(phi), 0, sind(phi);
    0, 1, 0;
    -sind(phi), 0, cosd(phi)];

R_z = [cosd(alpha), -sind(alpha), 0;
    sind(alpha), cosd(alpha), 0;
    0, 0, 1];

%--------------------------------
%  construct complete rotation matrix
%--------------------------------

A = R_x*R_y*R_z;

R_matrix = [A.^2,2*[A(1,2)*A(1,3), A(1,3)*A(1,1), A(1,1)*A(1,2); A(2,2)*A(2,3), A(2,3)*A(2,1), A(2,1)*A(2,2); A(3,2)*A(3,3), A(3,3)*A(3,1), A(3,1)*A(3,2)];
    A(2,1)*A(3,1), A(2,2)*A(3,2), A(2,3)*A(3,3), A(2,2)*A(3,3)+A(2,3)*A(3,2), A(2,1)*A(3,3)+A(2,3)*A(3,1), A(2,2)*A(3,1)+A(2,1)*A(3,2);
    A(3,1)*A(1,1), A(3,2)*A(1,2), A(3,3)*A(1,3), A(1,2)*A(3,3)+A(1,3)*A(3,2), A(1,3)*A(3,1)+A(1,1)*A(3,3), A(1,1)*A(3,2)+A(1,2)*A(3,1);
    A(1,1)*A(2,1), A(1,2)*A(2,2), A(1,3)*A(2,3), A(1,2)*A(2,3)+A(1,3)*A(2,2), A(1,3)*A(2,1)+A(1,1)*A(2,3), A(1,1)*A(2,2)+A(1,2)*A(2,1)];

%--------------------------------
%  rotate B into new frame
%--------------------------------

new_mtx = R_matrix*B*R_matrix';

    
    
