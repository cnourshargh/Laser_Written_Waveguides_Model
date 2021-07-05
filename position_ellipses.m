function [Txx_n,Tyy_n,Txy_n]=position_ellipses(Txx,Tyy,Txy,dx,positions)



g = length(Txx);

Txx_n = zeros(size(Txx));
Tyy_n = zeros(size(Tyy));
Txy_n = zeros(size(Txy));


for i = 1:2:length(positions)
    
    xposition = positions(i);
    yposition = positions(i+1);
    
    matrixshift_x = xposition/dx;
    
    matrixshift_y = yposition/dx;
    
    if xposition < 0 
        
        if yposition < 0
            Txx_n(1:g+matrixshift_y,1:g+matrixshift_x) = Txx_n(1:g+matrixshift_y,1:g+matrixshift_x) + Txx(-matrixshift_y+1:g,-matrixshift_x+1:g);
            Tyy_n(1:g+matrixshift_y,1:g+matrixshift_x) = Tyy_n(1:g+matrixshift_y,1:g+matrixshift_x) + Tyy(-matrixshift_y+1:g,-matrixshift_x+1:g);
            Txy_n(1:g+matrixshift_y,1:g+matrixshift_x) = Txy_n(1:g+matrixshift_y,1:g+matrixshift_x) + Txy(-matrixshift_y+1:g,-matrixshift_x+1:g);
        elseif yposition > 0
            Txx_n(matrixshift_y+1:g,1:g+matrixshift_x) = Txx_n(matrixshift_y+1:g,1:g+matrixshift_x) + Txx(1:g-matrixshift_y,-matrixshift_x+1:g);
            Tyy_n(matrixshift_y+1:g,1:g+matrixshift_x) = Tyy_n(matrixshift_y+1:g,1:g+matrixshift_x) + Tyy(1:g-matrixshift_y,-matrixshift_x+1:g);
            Txy_n(matrixshift_y+1:g,1:g+matrixshift_x) = Txy_n(matrixshift_y+1:g,1:g+matrixshift_x) + Txy(1:g-matrixshift_y,-matrixshift_x+1:g);
        else
            Txx_n(:,1:g+matrixshift_x) = Txx_n(:,1:g+matrixshift_x) + Txx(:,-matrixshift_x+1:g);
            Tyy_n(:,1:g+matrixshift_x) = Tyy_n(:,1:g+matrixshift_x) + Tyy(:,-matrixshift_x+1:g);
            Txy_n(:,1:g+matrixshift_x) = Txy_n(:,1:g+matrixshift_x) + Txy(:,-matrixshift_x+1:g);
        end
            
    elseif xposition > 0
        
        if yposition < 0
            Txx_n(1:g+matrixshift_y,matrixshift_x+1:g) = Txx_n(1:g+matrixshift_y,matrixshift_x+1:g) + Txx(-matrixshift_y+1:g,1:g-matrixshift_x);
            Tyy_n(1:g+matrixshift_y,matrixshift_x+1:g) = Tyy_n(1:g+matrixshift_y,matrixshift_x+1:g) + Tyy(-matrixshift_y+1:g,1:g-matrixshift_x);
            Txy_n(1:g+matrixshift_y,matrixshift_x+1:g) = Txy_n(1:g+matrixshift_y,matrixshift_x+1:g) + Txy(-matrixshift_y+1:g,1:g-matrixshift_x);
        elseif yposition > 0
            Txx_n(matrixshift_y+1:g,matrixshift_x+1:g) = Txx_n(matrixshift_y+1:g,matrixshift_x+1:g) + Txx(1:g-matrixshift_y,1:g-matrixshift_x);
            Tyy_n(matrixshift_y+1:g,matrixshift_x+1:g) = Tyy_n(matrixshift_y+1:g,matrixshift_x+1:g) + Tyy(1:g-matrixshift_y,1:g-matrixshift_x);
            Txy_n(matrixshift_y+1:g,matrixshift_x+1:g) = Txy_n(matrixshift_y+1:g,matrixshift_x+1:g) + Txy(1:g-matrixshift_y,1:g-matrixshift_x);
        else
            Txx_n(:,matrixshift_x+1:g) = Txx_n(:,matrixshift_x+1:g) + Txx(:,1:g-matrixshift_x);
            Tyy_n(:,matrixshift_x+1:g) = Tyy_n(:,matrixshift_x+1:g) + Tyy(:,1:g-matrixshift_x);
            Txy_n(:,matrixshift_x+1:g) = Txy_n(:,matrixshift_x+1:g) + Txy(:,1:g-matrixshift_x);
        end
    else
        
        if yposition < 0
            Txx_n(1:g+matrixshift_y,:) = Txx_n(1:g+matrixshift_y,:) + Txx(-matrixshift_y+1:g,:);
            Tyy_n(1:g+matrixshift_y,:) = Tyy_n(1:g+matrixshift_y,:) + Tyy(-matrixshift_y+1:g,:);
            Txy_n(1:g+matrixshift_y,:) = Txy_n(1:g+matrixshift_y,:) + Txy(-matrixshift_y+1:g,:);
        elseif yposition > 0
            Txx_n(matrixshift_y+1:g,:) = Txx_n(matrixshift_y+1:g,:) + Txx(1:g-matrixshift_y,:);
            Tyy_n(matrixshift_y+1:g,:) = Tyy_n(matrixshift_y+1:g,:) + Tyy(1:g-matrixshift_y,:);
            Txy_n(matrixshift_y+1:g,:) = Txy_n(matrixshift_y+1:g,:) + Txy(1:g-matrixshift_y,:);
        else
            Txx_n(:,:) = Txx_n(:,:) + Txx(:,:);
            Tyy_n(:,:) = Tyy_n(:,:) + Tyy(:,:);
            Txy_n(:,:) = Txy_n(:,:) + Txy(:,:);
        end
    end
    
    
end

