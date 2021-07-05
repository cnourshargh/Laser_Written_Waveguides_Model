function [] = meshplotter(gcoord,nodes)
%meshplotter.m produces a plot showing the positions of all the nodes and
%the outlines of all the finite elements.
%
%Inputs
%   gcoord = coordinates of all mesh nodes
%   nodes = matrix showing which nodes are in which finite element
%

sz = size(nodes);                       
n_elements = sz(1);                     %finds number of finite elements


figure;                                 %initialises figure

scatter(gcoord(:,1),gcoord(:,2),'x');   %plots all nodes
hold on;

for i = 1:n_elements                    %itterates over all finite elements
    plot(gcoord(nodes(i,[1,8,4,7,3,6,2,5,1]),1),gcoord(nodes(i,[1,8,4,7,3,6,2,5,1]),2),'b');
                                        %plots outline of element
    hold on;
end

xlabel('X coordinate (microns)');       %adds figure labels
ylabel('Y coordinate (microns)');
title('Mesh')
hold off;