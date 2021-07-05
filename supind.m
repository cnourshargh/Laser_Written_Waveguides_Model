function [Dnx,Dny,Dnz]=supind(Txx,Tyy, photoelastic, refractive_index,x,y,a,b,track_index,positions)
% Function to calculate change in refrative index profiles using
% photoelastic tensor
%
% INPUTS
% Txx ------------------ strain in x direction
% Tyy ------------------ strain in y direction
% photoelastic --------- photelastic matrix from material_selector
% refractive_index ----- refractive index matrix from material_selector
%
%
% OUTPUTS
% Dnx ------------------ change in refractive index for x polarisation
% Dny ------------------ change in refractive index for y polarisation
% Dnz ------------------ change in refractive index for z polarisation

if (nargin<9)
    track_index = 0;
    positions = [0,0];
else

end

Tzz = zeros(size(Txx));

Dnx = zeros(size(Txx));
Dny = zeros(size(Txx));
Dnz = zeros(size(Txx));

for i=1:length(Txx)
    for j=1:length(Tyy)
        
        Dnx(i,j)=-0.5*refractive_index(1,1)^3*(photoelastic(1,1)*Txx(i,j)+photoelastic(1,2)*Tyy(i,j)+photoelastic(1,3)*Tzz(i,j));   
        Dny(i,j)=-0.5*refractive_index(2,2)^3*(photoelastic(2,1)*Txx(i,j)+photoelastic(2,2)*Tyy(i,j)+photoelastic(2,3)*Tzz(i,j));
        Dnz(i,j)=-0.5*refractive_index(3,3)^3*(photoelastic(3,1)*Txx(i,j)+photoelastic(3,2)*Tyy(i,j)+photoelastic(3,3)*Tzz(i,j));
        
        for k = 1:length(positions)/2
            if ((x(i,j)-positions(2*k-1))/a)^2+((y(i,j)-positions(2*k))/b)^2 < 0.9999                   %if coordinate is within damage elipse
                Dnx(i,j)=-0.00+track_index;                         %set refractive index change to pure imaginary
                Dny(i,j)=-0.00+track_index;
                Dnz(i,j)=-0+track_index;
            end
        end
    end
end
        

