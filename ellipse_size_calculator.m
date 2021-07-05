function [a,b] = ellipse_size_calculator(lambda,NA)
% Function to calculate ellipse geometry using Abbe Diffraction Limit

% INPUTS
% lambda ---- Peak wavelength (in microns) of machining laser
% NA -------- Numerical aperture of machining objective

% OUPUTS
% a --------- Ellipse radius in x
% b --------- Ellipse radius in y

    a = round(lambda/(2*NA),2,'significant');              % Lateral resolution gives x radius
    b = round(2*lambda/(NA^2),2,'significant');            % Axial resolution gives y radius

end