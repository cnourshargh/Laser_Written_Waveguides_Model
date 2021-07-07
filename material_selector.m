function [stiffness, photoelastic, refractive_index] = material_selector(opt,rotation)
%material_selector.m creates the stiffness, photoelastic, and refractive
%indx matrices required for the calculation of the refractive index change.
%Built in are some salient cuts of lithium niobate, and options for user
%defined orientations of lithium niobate and sapphire. More materials can
%be added by putting the (unrotated) matrices in the 'Defining material
%properties' section and adding an option for them to be called in the
%'Select Material' section.
%
%Inputs
%   opt = option from material library
%   rotation = vector describing matrix rotation, of form [theta, phi, alpha].
%               For any arbitrary cut options, this must be defined.        
%   theta = rotation of material matrices about the x axis
%   phi = rotation of material matrices about the y axis
%   alpha = rotation of material matrices about the z axis
%
%Outputs
%   Stiffness = 3x3 matrix describing material stifness under plane strain
%                   conditions
%   photoelastic = 3x3 matrix of pgotoelastic coefficients relevant to the
%                   material orientation being considered
%   refracive index = refractive index matrix, alligned with the
%                       orientation of the material
%
%Preprogrammed Options
%   1 - Lithium Niobate, x cut, horizontal z, y out of plane
%   2 - Lithium Niobate, x cut, vertical z, y out of plane
%   3 - Lithium Niobate arbitrary cut
%   4 - Sapphire arbitrary cut

%--------------------------------
%  Defining material properties
%--------------------------------

LiNbO3_stiffness = [0.2030  0.0573  0.0752  0.0085  0  0;
                    0.0573  0.2030  0.0752  -0.0085  0  0;
                    0.0752  0.0752  0.2424  0  0  0;
                    0.0085  -0.0085  0  0.0595  0  0;
                    0  0  0  0  0.0595  0.0085;
                    0  0  0  0  0.0085  0.0729]; 
                %3D stiffness matrix, coefficients given in Newtons per
                %square micron (10^12 Pa)
                %taken from https://link.springer.com/article/10.1007/BF00614817
                
                
LiNbO3_photoelastic = [-0.026  0.090  0.133  -0.075  0  0;
                        0.090  -0.026  0.133  0.075  0  0;
                        0.179  0.179  0.071  0  0  0;
                        -0.151  0.151  0  0.146  0  0;
                        0  0  0  0  0.146  -0.151;
                        0  0  0  0  -0.075  -0.058];
                %3D photoelastic matrix at 633nm, taken from Newnham
                
LiNbO3_index = [2.2580, 0, 0;
                0, 2.2580, 0;
                0, 0, 2.1778];
                %refractive inex matrix, taken from Zelmon et al. 1997
                
Sapphire_stiffness = [0.486  0.098  0.063  0.037  0  0;
                      0.098  0.486  0.063  -0.037  0  0;
                      0.063  0.063  0.496  0  0  0;
                      0.037  -0.037  0  0.186  0  0;
                      0  0  0  0  0.186  0.037;
                      0  0  0  0  0.037  0.194];
                  %stiffness matrix for sapphire, given in given in Newtons per
                  %square micron (10^12 Pa)
                  % taken from  https://doi.org/doi:10.25335/M5782R
                
Sapphire_photoelastic = [-0.23  -0.03  0.02  0.00  0  0;
                         -0.03  -0.23  0.02  -0.00  0  0;
                         -0.04  -0.04  -0.20  0  0  0;
                         0.01  -0.01  0  -0.10  0  0;
                         0  0  0  0  -0.10  0.01;
                         0  0  0  0  0.00  0.10];   
                 %photoelastic matrix for sapphire taken from newnham 

Sapphire_index = [1.7601  0  0;
                  0  1.7601  0;
                  0  0  1.7522];
              %refractive index matrix for sapphire, at 800nm
              %taken from https://refractiveindex.info/?shelf=main&book=Al2O3&page=Malitson-o

              
%--------------------------------
%  Selecting Material
%--------------------------------

if opt == 1                 %given option 1 - Lithium Niobate, x cut, horizontal z, y out of plane             
    
    rotation = [90,90,0];   %rotation is predefined
    
    stiffness_int = tensor_rotator(LiNbO3_stiffness,rotation);
                            %6x6 stiffness matrix is rotated into lab frame
    
    stiffness = stiffness_int([1,2,6],[1,2,6]);
                            %relevant values for plane strain assumption
                            %are identified

    phot_int = tensor_rotator(LiNbO3_photoelastic,rotation);
                            %6x6 photoelastic matrix is rotated into lab frame
    photoelastic = phot_int([1,2,6],[1,2,6]);
                            %relevant values for plane strain assumption
                            %are identified
    
    refractive_index = matrix_rotator(LiNbO3_index,rotation);
                            %3x3 refracive index matrix is rotated into lab
                            %frame
    
elseif opt==2
    
    rotation = [90,0,0];
    
    stiffness_int = tensor_rotator(LiNbO3_stiffness,rotation);
    
    stiffness = stiffness_int([1,2,6],[1,2,6]);

    phot_int = tensor_rotator(LiNbO3_photoelastic,rotation);
    photoelastic = phot_int([1,2,6],[1,2,6]);
    
    refractive_index = matrix_rotator(LiNbO3_index,rotation);
    
    
elseif opt==3
    
    stiffness_int = tensor_rotator(LiNbO3_stiffness,rotation);
    
    stiffness = stiffness_int([1,2,6],[1,2,6]);

    phot_int = tensor_rotator(LiNbO3_photoelastic,rotation);
    photoelastic = phot_int([1,2,6],[1,2,6]);
    
    refractive_index = matrix_rotator(LiNbO3_index,rotation);
    
    
elseif opt==4
    
    stiffness_int = tensor_rotator(Sapphire_stiffness,rotation);
    
    stiffness = stiffness_int([1,2,6],[1,2,6]);

    phot_int = tensor_rotator(Sapphire_photoelastic,rotation);
    photoelastic = phot_int([1,2,6],[1,2,6]);
    
    refractive_index = matrix_rotator(Sapphire_index,rotation);
    
    
end



end

                    