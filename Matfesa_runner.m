function [] = Matfesa_runner()
%Executes subroutines of Matfesa code to use finite element method to
%caluculate stress induced refracive index profile, and then finite
%difference method to estimate to intensity profiles of modes guided by
%laser written waveguides
%
%Inputs:
%   a,b = track X- and Y- radii respectively
%   ellipticity = n in equation (x/a)^n + y(/b)^n = 1 used to define ellipse
%   linear_expansion = linear expansion coefficient applied to strain (track) boundary
%   quadratic_expansion = quadratic expansion coefficient applied to strain (track) boundary               
%   positions = vector describing positions of tracks in two dimensions
%               [x1,y1,x2,y2,...]
%   opt = Used to select material parameters, see material_selector.m for
%       more
%   rotation = used to define material orientation when not specified in
%               material_selector.m
%   track_index = refractive index change applied within the focal volume
%   lambda_illumination = wavelength of light considered when calulating
%                       guided modes
%   lambda_machining, NA_machining = wavelength and numerical apperture of
%   laser machining system, used to calculate ellipse size in Abbe limitted
%   case


%------------------------------------
%define Inputs
%------------------------------------

linear_expansion = 0.01;                    %linear expansion parameter
quadratic_expansion = 0;                    %quadratic expansion parameter

lambda_machining = .8;                      %Peak wavelength of machining laser
NA_machining = 0.5;                         %Numerical aperture of machining system
%[a,b] = ellipse_size_calculator(lambda_machining,NA_machining);
%calculates track dimensions from Abbe limit

a=1;
b=5;
positions = [0,0];  
%positions = [8,0,-8,0];
%positions = [5,0,-5,0,0,10,0,-10];
%positions = [-8,0,8,0,-6,8,6,8,-6,-8,6,-8,4,10,-4,10,-4,-10,4,-10,2,12,2,-12,-2,-12,-2,12,0,12,0,-12];

ellipticity = 2;                            %index of (x/a)^n + y(/b)^n = 1 in ellpse size
track_index = 0.00;                         %absorption coefficient (1/m) of damage tracks   50000

opt = 1;                                    %defines material and orientation
                                            %see material_selector.m for details
rotation = [90,90,0];                       %crystal orientation 
                                            %only used for unsepcified orientation
                                            %option
lambda_illumination = .78;                  %illumination wavelength in microns for testing
    
%-----------------------------
%calculates position of erroneous element produced by matfesa
%-----------------------------
[bad_element_boundary_y]=matfesa_simple_mesh(linear_expansion,quadratic_expansion,linear_expansion,quadratic_expansion,a,b,ellipticity,opt,rotation);

%-----------------------------
%calculates index profile for given track size and position
%-----------------------------
[nx,ny,nz,x,y,RI,dx]=matfesa_adaptive_mesh(linear_expansion,quadratic_expansion,linear_expansion,quadratic_expansion,a,b,bad_element_boundary_y,ellipticity,positions,opt,rotation,track_index);

%-----------------------------
%remove lingering peaks from real part of index profile with median
%filter
%-----------------------------
ny = medfilt2(real(ny))+1i*imag(ny);
nx = medfilt2(real(nx))+1i*imag(nx);            
nz = medfilt2(real(nz))+1i*imag(nz);

%-----------------------------
%remove outer rows and collumns from RI profiles and coordinate meshs
%as these are corrupted by the medfilt
%-----------------------------
nx = nx(2:end-1,2:end-1);                       
ny = ny(2:end-1,2:end-1);
nz = nz(2:end-1,2:end-1);
xp = x(2:end-1,2:end-1);
yp = y(2:end-1,2:end-1);

%-----------------------------
%plot real part of index profiles
%-----------------------------
figure;
surf(xp,yp,real(ny),'linestyle','none');
colorbar;
view(2);
title('Refractive Index Profile in Y');
xlabel('X coordinate (microns)');
ylabel('Y coordinate (microns)');

figure;
surf(xp,yp,real(nx),'linestyle','none');
colorbar;
view(2);
xlabel('X coordinate (microns)');
ylabel('Y coordinate (microns)');
title('Refractive Index Profile in X')

%-----------------------------
% finds ordinary and extraordinary refractive indices from refractive index
% matrix
%-----------------------------
seed_max = max(max(RI));
seed_min = min(RI(RI>0));    

%-----------------------------
% finds modes guided bthe polarisation with the highest unstrained index
% if matrial is naturally birefringent, also calculates modes guided by lowest index 
%-----------------------------
[Ix_high,Iy_high,neff_high]=find_wg_modes(RI,lambda_illumination,nx,ny,nz,dx,x,y,seed_max);

if seed_max ~= seed_min

    [Ix_low,Iy_low,neff_low]=find_wg_modes(RI,lambda_illumination,nx,ny,nz,dx,x,y,seed_min);

end


