function [gcoord,nodes]=mesh1(ao,bo,plot_size,mesh_res,bad_element_coord_boundaries,ellipticity)
%mesh1.m produces the finite element mesh for both the simple and adaptive
%mesh programs. It also generates the connectivities of all the nodes, used
%when calculating the porogation of stresses throughout the material. 
%
%Inputs
%   ao,bo = Boundary ellipse X- and Y- radii
%   plot_size = side length of region over which finite element model is
%               needed. Here this length is doubled to ensure the finite
%               element mesh is fine in the area around the tracks.
%   mesh_res = resolution metric for mesh row spacings, used as exponent coefficient 
%               in radial exponential 
%   bad_element_coord_boundaries = used in matfesa_adaptive_mesh.m.
%                                   Coordinates on the strain boundary
%                                   between which a higher resolution is
%                                   required.
%   ellipticity = n in equation (x/a)^n + y(/b)^n = 1 used to define strain
%               boundary
%Outputs
%   gcoord = collumn vector of x and  coordinates of all mesh nodes
%   nodes = matrix showing connections of each of the nodes in the finite
%           element mesh. Each row represents a finite element, made of
%           nine nodes

%-----------------------------------------
% Calculate number of element rows required for area
%-----------------------------------------

n_rows = 2;                         %number of rows of nodes required, 
                                    %calculated below

while ao+1*exp((mesh_res*n_rows))<2*plot_size
                                    %row width increases exponetiall with
                                    %row number until a is twice the plot
                                    %size
    
    n_rows = n_rows + 2;            %counts the number of rows required
                                    
end

%-----------------------------------------
% Find node coordinates
%-----------------------------------------

for lo=0:n_rows                     %itterates over all rows
    a = ao*exp((mesh_res*lo));      %calculates row X- radius
    b = bo+(a-ao);                  %calculates row Y- radius
    gco = fun1(a,b,bad_element_coord_boundaries,bo,ellipticity);
                                    %calculates coordinates of all nodes
                                    %along row
                                    
    iin = 1+lo*length(gco);         %finds start and end index for 
    ifi = iin+(length(gco)-1);      %coordinates in gcoord
    
    gcoord(iin:ifi,1:2) = gco;      %populates gcoord with node coordinates
end


%-----------------------------------------
% Define nodal connectivities
%-----------------------------------------

no = length(gco);                       %number of nodes per row
ftot = lo/2;                            %total number of rows of elements
nodes = zeros(((no-1)/2)*ftot,9);       %defines nodes matrix size and populates it with zeros
w = 0;                                  %counter to be used when iterating rows

for r = 1:ftot                          %itterating over all element rows
    iel = (r - 1)*(no - 1)/2 + 1;       %Finda index of first node in element row      
    for s = iel:iel + (no - 3)/2        %itterates over all elements in element row
    if s == iel                         %if first element in row
        
        i = w * no + 1;                 %identifies index of all nodes in element
        j = i + 2*no;
        k = j + 2;
        l = i + 2;
        m = i + length(gco);
        n = j + 1;
        o = l + length(gco);
        p = i + 1;
        q = m + 1;
        nodes(s,:) = [i j k l m n o p q];   
                                        %saves indexes in nodes matrix,
    
    else
        i = i + 2;                      %identifies index of all nodes in element
        j = i + 2*no;
        k = j + 2;
        l = i + 2;
        m = i + length(gco);
        n = j + 1;
        o = l + length(gco);
        p = i + 1;
        q = m + 1;
        nodes(s,:) = [i j k l m n o p q];
                                        %saves indexes in nodes matrix,
    end
    end
    
    w = w + 2;                          %itterates counter
    
end