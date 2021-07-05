function [Txx_n,Tyy_n,Txy_n,x_n,y_n] = plot_full_strain_map(Txx,Tyy,Txy,dx)
    
    g = length(Txx);
    h = (g-1)*dx;
    
    Txx_n = zeros(2*g-1,2*g-1);
    Tyy_n = zeros(2*g-1,2*g-1);
    Txy_n = zeros(2*g-1,2*g-1);
    
    %bottom left quadrant
    Txx_n(1:g,1:g) = Txx(g:-1:1,g:-1:1);
    Tyy_n(1:g,1:g) = Tyy(g:-1:1,g:-1:1);
    Txy_n(1:g,1:g) = Txy(g:-1:1,g:-1:1);
    
    %bottom right quadrant
    Txx_n(1:g,g:end) = Txx(g:-1:1,1:g);
    Tyy_n(1:g,g:end) = Tyy(g:-1:1,1:g);
    Txy_n(1:g,g:end) = Txy(g:-1:1,1:g);
    
    %top left quadrant
    Txx_n(g:end,1:g) = Txx(1:g,g:-1:1);
    Tyy_n(g:end,1:g) = Tyy(1:g,g:-1:1);
    Txy_n(g:end,1:g) = Txy(1:g,g:-1:1);
    
    %bottom right quadrant
    Txx_n(g:end,g:end) = Txx(1:g,1:g);
    Tyy_n(g:end,g:end) = Tyy(1:g,1:g);
    Txy_n(g:end,g:end) = Txy(1:g,1:g);
    
    [x_n,y_n] = meshgrid(-h:dx:h,-h:dx:h);
        
end