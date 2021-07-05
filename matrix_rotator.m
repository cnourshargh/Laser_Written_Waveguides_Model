function new_mtx = matrix_rotator(B,rotation)

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

new_mtx = A*B*A';