function new_mtx = matrix_rotator(B,rotation)
%matrix_rotator.m is used to rotate the refractive index matrix in 
%material_selector.m. It constructs 3 rotation matrices, relating to tjhe
%three axes about which rotations are defined and multiplies all three
%together, and then rotates refractie index matrix into the new frame.
%
%Inputs
%   B = 3x3 matrix to be rotated
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

%--------------------------------
%  Rotate B into the frame defined by the rotation
%--------------------------------

new_mtx = A*B*A';