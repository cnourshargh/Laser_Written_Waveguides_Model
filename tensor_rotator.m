function new_mtx = tensor_rotator(B,rotation)

theta = rotation(1);
phi = rotation(2);
alpha = rotation(3);

R_x = [1, 0, 0;
    0, cosd(theta), -sind(theta);
    0, -sind(theta), cosd(theta)];

R_y = [cosd(phi), 0, sind(phi);
    0, 1, 0;
    -sind(phi), 0, cosd(phi)];

R_z = [cosd(alpha), -sind(alpha), 0;
    sind(alpha), cosd(alpha), 0;
    0, 0, 1];

A = R_x*R_y*R_z;

R_matrix = [A.^2,2*[A(1,2)*A(1,3), A(1,3)*A(1,1), A(1,1)*A(1,2); A(2,2)*A(2,3), A(2,3)*A(2,1), A(2,1)*A(2,2); A(3,2)*A(3,3), A(3,3)*A(3,1), A(3,1)*A(3,2)];
    A(2,1)*A(3,1), A(2,2)*A(3,2), A(2,3)*A(3,3), A(2,2)*A(3,3)+A(2,3)*A(3,2), A(2,1)*A(3,3)+A(2,3)*A(3,1), A(2,2)*A(3,1)+A(2,1)*A(3,2);
    A(3,1)*A(1,1), A(3,2)*A(1,2), A(3,3)*A(1,3), A(1,2)*A(3,3)+A(1,3)*A(3,2), A(1,3)*A(3,1)+A(1,1)*A(3,3), A(1,1)*A(3,2)+A(1,2)*A(3,1);
    A(1,1)*A(2,1), A(1,2)*A(2,2), A(1,3)*A(2,3), A(1,2)*A(2,3)+A(1,3)*A(2,2), A(1,3)*A(2,1)+A(1,1)*A(2,3), A(1,1)*A(2,2)+A(1,2)*A(2,1)];

new_mtx = R_matrix*B*R_matrix';

    
    
