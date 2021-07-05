function [dnx_max,dny_max,dnx_min,dny_min,epsilon] = index_vs_expansion()

n_datapoints = 3;

dnx_max = zeros([1,n_datapoints]);
dny_max = zeros([1,n_datapoints]);
dnx_min = zeros([1,n_datapoints]);
dny_min = zeros([1,n_datapoints]);

% dnx_max_true = zeros([1,n_datapoints]);
% dnx_min_true = zeros([1,n_datapoints]);
% dny_max_true = zeros([1,n_datapoints]);
% dny_min_true = zeros([1,n_datapoints]);

epsilon = logspace(-5,-2,n_datapoints);


lambda_illumination = .78;  

quadratic_expansion = 0;  

positions = [0,0];

track_absorption = 0;

a=1;
b=1;

opt = 1;                            %defines material and orientation
                                    %see material_selector.m for details
rotation = [];                      %crystal orientation 
                                    %only used for unsepcified orientation
                                    %option
ellipticity =2;                     %index of (x/a)^n + y(/b)^n = 1 in ellpse size

for i = 1:n_datapoints
    
    %-----------------------------
    %calculates position of erroneous element produced by matfesa
    %-----------------------------
    [bad_element_boundary_y]=matfesa_simple_mesh(epsilon(i),quadratic_expansion,epsilon(i),quadratic_expansion,a,b,ellipticity,opt,rotation);
    %-----------------------------
    %calculates index profile for given track size and position
    %-----------------------------
    [nx,ny,nz,x,y,RI,dx,Dnx,Dny,Dnz,PE]=matfesa_adaptive_mesh(epsilon(i),quadratic_expansion,epsilon(i),quadratic_expansion,a,b,bad_element_boundary_y,ellipticity,positions,opt,rotation,track_absorption);
    %-----------------------------
    %remove lingering peaks from real part of index profile with median
    %filter
    %-----------------------------
    Dnx = medfilt2(real(Dnx))+1i*imag(Dnx);            
    Dny = medfilt2(real(Dny))+1i*imag(Dny);
    Dnz = medfilt2(real(Dnz))+1i*imag(Dnz);
    %-----------------------------
    %remove outer rows and collumns from RI profiles and coordinate meshs
    %as these are corrupted by the medfilt
    %-----------------------------
    Dnx = Dnx(2:end-1,2:end-1);                       
    Dny = Dny(2:end-1,2:end-1);
    Dnz = Dnz(2:end-1,2:end-1);
    xp = x(2:end-1,2:end-1);
    yp = y(2:end-1,2:end-1);
    
    dnx_max(i) = max(max(Dnx(:,:)));
    dnx_min(i) = min(min(Dnx(:,:)));
    dny_max(i) = max(max(Dny(:,:)));
    dny_min(i) = min(min(Dny(:,:)));

    
%     dnx_max_true(i) = (RI(1,1)^(-2)+epsilon(i)*PE(1,2)*b)^(-0.5)-RI(1,1);
%     dnx_min_true(i) = (RI(1,1)^(-2)+epsilon(i)*PE(1,1)*a)^(-0.5)-RI(1,1);
%     dny_max_true(i) = (RI(2,2)^(-2)+epsilon(i)*PE(2,1)*a)^(-0.5)-RI(2,2);
%     dny_min_true(i) = (RI(2,2)^(-2)+epsilon(i)*PE(2,2)*b)^(-0.5)-RI(2,2);
end

figure;
plot(epsilon,dnx_max)
xlabel('Boundary Strain');
ylabel('Maximum Change in X- Index');
title('Maximum X Index vs Expansion Parameter')

figure;
plot(epsilon,dnx_min)
xlabel('Expansion Parameter');
ylabel('Minimum Value of dNx');
title('Minimum X Index vs Expansion Parameter')

figure;
plot(epsilon,dny_max)
xlabel('Expansion Parameter');
ylabel('Maximum Value of dNy');
title('Maximum Y Index vs Expansion Parameter')

figure;
plot(epsilon,dny_min)
xlabel('Expans ion Parameter');
ylabel('Minimum Value of dNy');
title('Minimum Y Index vs Expansion Parameter')