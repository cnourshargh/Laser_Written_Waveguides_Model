function [F_Ix,F_Iy,F_Ix_var,F_Iy_var,F_I_var] = fourier_spectral_varience(Ix,Iy,neff,seed,x,y,n_modes)

%-----------------------------
%Initailise matrices for looping
%-----------------------------

F_Ix = zeros(size(Ix));
F_Iy = zeros(size(Ix));
F_Ix_normalised = zeros(size(Ix));
F_Iy_normalised = zeros(size(Ix));

F_Ix_mean = zeros([2,n_modes]);
F_Iy_mean = zeros([2,n_modes]);

F_Ix_var = zeros([2,n_modes]);
F_Iy_var = zeros([2,n_modes]);

F_I_var = zeros([2,n_modes]);

for i = 1:n_modes
    
    %-----------------------------
    %calculate normalised magnitude fourier spectrum of mode profiles
    %-----------------------------
    
    F_Ix(:,:,i) = abs(fftshift(fft2(Ix(:,:,i))));
    F_Iy(:,:,i) = abs(fftshift(fft2(Iy(:,:,i))));
    
    F_Ix_normalised(:,:,i) = (1/max(max(F_Ix(:,:,i)))) * F_Ix(:,:,i);
    F_Iy_normalised(:,:,i) = (1/max(max(F_Iy(:,:,i)))) * F_Iy(:,:,i);
   
    %-----------------------------
    %calculate centroids of fourier spectra
    %-----------------------------
    
    F_Ix_mean(1,i) = (1/sum(F_Ix_normalised(:,:,i),'all'))*sum(x.*F_Ix_normalised(:,:,i),'all');
    F_Ix_mean(2,i) = (1/sum(F_Ix_normalised(:,:,i),'all'))*sum(y.*F_Ix_normalised(:,:,i),'all');
    
    F_Iy_mean(1,i) = (1/sum(F_Iy_normalised(:,:,i),'all'))*sum(x.*F_Iy_normalised(:,:,i),'all');
    F_Iy_mean(2,i) = (1/sum(F_Iy_normalised(:,:,i),'all'))*sum(y.*F_Iy_normalised(:,:,i),'all');
    
    %-----------------------------
    %calculate spatial variance of fourier spectra
    %-----------------------------
    
    F_Ix_var(1,i) =  (1/sum(F_Ix_normalised(:,:,i),'all'))*sum((x.^2).*F_Ix_normalised(:,:,i),'all')-F_Ix_mean(1,i).^2;
    F_Ix_var(2,i) =  (1/sum(F_Ix_normalised(:,:,i),'all'))*sum((y.^2).*F_Ix_normalised(:,:,i),'all')-F_Ix_mean(2,i).^2;
    
    F_Iy_var(1,i) =  (1/sum(F_Iy_normalised(:,:,i),'all'))*sum((x.^2).*F_Iy_normalised(:,:,i),'all')-F_Iy_mean(1,i).^2;
    F_Iy_var(2,i) =  (1/sum(F_Iy_normalised(:,:,i),'all'))*sum((y.^2).*F_Iy_normalised(:,:,i),'all')-F_Iy_mean(2,i).^2;
    
    F_I_var(1,i) = F_Ix_var(1,i) + F_Ix_var(2,i);
    F_I_var(2,i) = F_Iy_var(1,i) + F_Iy_var(2,i);
    
    if neff(i) > seed
        guided(i) = 1;
    else
        guided(i) = 0;
end
end

for j = 1:n_modes
if guided(j) == 1
figure;
subplot(2,2,1);
surf(x,y,Iy(:,:,j),'linestyle','none');
view(2);
colorbar;
xlabel('X coordinate (microns)');
ylabel('Y coordinate (microns)');
title('Mode Profile')

subplot(2,2,2);
surf(x,y,F_Iy(:,:,j),'linestyle','none');
view(2);
colorbar;
xlabel('X coordinate (microns)');
ylabel('Y coordinate (microns)');
title('Fourier Spectrum')

subplot(2,2,3);
surf(x,y,Ix(:,:,j),'linestyle','none');
view(2);
colorbar;
xlabel('X coordinate (microns)');
ylabel('Y coordinate (microns)');
title('Mode Profile')

subplot(2,2,4);
surf(x,y,F_Ix(:,:,j),'linestyle','none');
view(2);
colorbar;
xlabel('X coordinate (microns)');
ylabel('Y coordinate (microns)');
title('Fourier Spectrum')
end
end

[m,j] = min(min(abs(64-F_I_var)));


figure;
subplot(2,2,1);
surf(x,y,Iy(:,:,j),'linestyle','none');
view(2);
colorbar;
xlabel('X coordinate (microns)');
ylabel('Y coordinate (microns)');
title('Mode Profile')

subplot(2,2,2);
surf(x,y,F_Iy(:,:,j),'linestyle','none');
view(2);
colorbar;
xlabel('X coordinate (microns)');
ylabel('Y coordinate (microns)');
title('Fourier Spectrum')

subplot(2,2,3);
surf(x,y,Ix(:,:,j),'linestyle','none');
view(2);
colorbar;
xlabel('X coordinate (microns)');
ylabel('Y coordinate (microns)');
title('Mode Profile')

subplot(2,2,4);
surf(x,y,F_Ix(:,:,j),'linestyle','none');
view(2);
colorbar;
xlabel('X coordinate (microns)');
ylabel('Y coordinate (microns)');
title('Fourier Spectrum')