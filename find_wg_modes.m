function [Ix,Iy,neff] = find_wg_modes(RI,lambda_illumination,nx,ny,nz,dx,x,y,seed)
%takes index profiles and finds fundmental guided modes 
%if material is originally birefrigent, runs modeslver twice to find modes
%guided around each ambient index
%
% RI ------------------- matrix of unstrained refractive indeces
% Lambda_illumniation -- Wavelength of light in microns propogating through sturcture
% nx, ny, nz ----------- refractive index profiles 
% dx ------------------- mesh size in microns

%-----------------------------
%run modesolver with seedmax and square field for intensity profile
%-----------------------------
n_modes = 40;
[hx,hy,neff] = wgmodes(lambda_illumination,seed,n_modes,dx,0.2,nx.^2,ny.^2,nz.^2,'0000');

Ix = abs(hx).^2;
Iy = abs(hy).^2;
%-----------------------------
%Initailise matrices for looping
%-----------------------------
x = x(2:end,2:end)-dx/2;
y = y(2:end,2:end)-dx/2;

Ix_normalised = zeros(size(Ix));
Iy_normalised = zeros(size(Iy));

Ix_mean = zeros([2,size(neff)]);
Iy_mean = zeros([2,size(neff)]);

Ix_var = zeros([2,size(neff)]);
Iy_var = zeros([2,size(neff)]);

I_var =  zeros([2,size(neff)]);

guided = zeros([size(neff),1]);


for i = 1:size(neff)
        
    %-----------------------------
    %calculate normalised intensity profile spatial variances
    %-----------------------------
    
    Ix_normalised(:,:,i) = (1/max(max(Ix(:,:,i)))) * Ix(:,:,i);
    Iy_normalised(:,:,i) = (1/max(max(Iy(:,:,i)))) * Iy(:,:,i);
    
    %-----------------------------
    %calculate centre of intensity profile as spatial mean
    %-----------------------------
    
    Ix_mean(1,i) = (1/sum(Ix_normalised(:,:,i),'all'))*sum(x.*Ix_normalised(:,:,i),'all');
    Ix_mean(2,i) = (1/sum(Ix_normalised(:,:,i),'all'))*sum(y.*Ix_normalised(:,:,i),'all');
    
    Iy_mean(1,i) = (1/sum(Iy_normalised(:,:,i),'all'))*sum(x.*Iy_normalised(:,:,i),'all');
    Iy_mean(2,i) = (1/sum(Iy_normalised(:,:,i),'all'))*sum(y.*Iy_normalised(:,:,i),'all');
    
    %-----------------------------
    %calculate variance of intensity profile 
    %-----------------------------
    
    Ix_var(1,i) =  (1/sum(Ix_normalised(:,:,i),'all'))*sum((x.^2).*Ix_normalised(:,:,i),'all')-Ix_mean(1,i).^2;
    Ix_var(2,i) =  (1/sum(Ix_normalised(:,:,i),'all'))*sum((y.^2).*Ix_normalised(:,:,i),'all')-Ix_mean(2,i).^2;
    
    Iy_var(1,i) =  (1/sum(Iy_normalised(:,:,i),'all'))*sum((x.^2).*Iy_normalised(:,:,i),'all')-Iy_mean(1,i).^2;
    Iy_var(2,i) =  (1/sum(Iy_normalised(:,:,i),'all'))*sum((y.^2).*Iy_normalised(:,:,i),'all')-Iy_mean(2,i).^2;
    
    I_var(1,i) = Ix_var(1,i) + Ix_var(2,i);
    I_var(2,i) = Iy_var(1,i) + Iy_var(2,i);
    
    %-----------------------------
    %determine ifmode is guided
    %-----------------------------
    
    if neff(i) > seed
        guided(i) = 1;
    else
        guided(i) = 0;
    end 
end
[m,j] = min(min(I_var));
for j = 1:n_modes
if guided(j) == 1 && min(I_var(:,j))<200
    
    caption = "Effective index = "+neff(j); 
    dim = [0.5, 0.2, 0.1, 0.1];

    figure;
    surf(x,y,Iy(:,:,j),'linestyle','none');
    view(2);
    colorbar;
    xlabel('X coordinate (microns)');
    ylabel('Y coordinate (microns)');
    title('Mode Intensity profile for X- Polarisation State')
    %annotation('textbox',dim,'String',caption);

    figure;
    surf(x,y,Ix(:,:,j),'linestyle','none');
    view(2);
    colorbar;
    xlabel('X coordinate (microns)');
    ylabel('Y coordinate (microns)');
    title('Mode Intensity profile for Y- Polarisation State')
    %annotation('textbox',dim,'String',caption);
end
end
figure;
scatter(real(neff(:)),I_var(1,:))
xlabel('Mode effective index');
ylabel('Mode Spatial Variance');
title('Spatial variacne vs Effective index')
% 
figure;
scatter(real(neff(:)),I_var(2,:))
xlabel('Mode effective index');
ylabel('Mode Spatial Variance');
title('Spatial variacne vs Effective index')
end