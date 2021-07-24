% clc
clear
close all
% z = im2double(imread('5_barbara.png'));
% z=im2double(imread('barbara.tif'));
% z = im2double(imread('6_boat.png'));
% z = im2double(imread('Lena.png'));
z = im2double(imread('lena256.tif'));
% z = im2double(imread('Cameraman.png'));
% z = im2double(imread('2_house.png'));
% z = im2double(imread('house.tif'));
% z=im2double(imread('barbara.png'));
% z=im2double(imread('peppers.png'));
% % z = im2double(imread('boats.tif'));
% z=im2double(imread('peppers.tif'));
% z=im2double(imread('cameraman.tif'));
%  z=im2double(imread('Monarch_full.tif'));
%   s1=0; for a1=-7:7; s1=s1+1; s2=0; for a2=-7:7; s2=s2+1; h(s1,s2)=1/(a1^2+a2^2+1); end, end;  h=h./sum(h(:));

% h = fspecial('gaussian',[9 9],1.2);%%%模糊核设置
h = fspecial('average',[9 9]);%%%模糊核设置
% h = fspecial('gaussian',[9 9]);%%%模糊核设置
% h = fspecial('gaussian',9, 1.1);
% z=im2double(imread('2_house.png'));
H=h;
x_orig=z;
noise_level = 2/255;%%

b = imfilter(z,h,'circular')+noise_level*randn(size(z));
b(b>1) = 1;b(b<0) = 0;
x0=b;%%初始化
figure,imshow(x0,[])

Hty = cconv2_by_fft2(b,h,0,[]);

gradf=@(x) HtH(x)-Hty;




searchRad = 10;             % Search window radius in NLM
patchRad = 3;               % Patch radius in NLM
hh=5/255;                 % Gaussian parameter in (SIGMA=20//NLM210) lena(6.5)



tol = 1e-5;                 % Tolerance for termination of algorithm
maxiters = 100;             % Max. no. of iterations

% delta = [0.3,0.6,0.9,1.0,1.2,1.5,1.7,1.9,2.0,2.1];
% delta = [0.3,0.5,0.7,0.9,1.2];
% delta = [0.6,0.9,1.2,1.6,1.9];
delta =0.9;
% % delta = [0.3,0.5,0.9,1.9];
%% Main algorithm
x0=b;
% x0 = initInpainting(b,P==0,5);  % Initial point to start the iterations
% W= @(x) staticBM3D(0, x, sigma, 'np', 0 ,x0);
figure; imshow(x0); title('Initialization');
psnr(x0,x_orig)


p2=psnr(x0,x_orig)
currentFile2 = sprintf('x0psnr_impaintingnlmx10%d',p2)
folder_name2 = ['C:\Users\47058\Desktop\2020pnp\实验数据\' currentFile2  '.png']
imwrite(x0, folder_name2);


Spsnrs=zeros(length(delta),maxiters);
Dpsnrs=zeros(length(delta),maxiters);
Serrors=zeros(length(delta),maxiters);
Derrors=zeros(length(delta),maxiters);
% sigma=10;
% sigma=hh;%%debluring
%dynamic


W1 = @(x) GNLM(x,x0,patchRad,searchRad,hh); 


for i=1:1:length(delta)
[x_hat1,converged,iters,errors,psnrs] = ...
        pnpISTA4(x0,gradf,W1,delta(i),x_orig,tol,10);  % Run algorithm
% x0=x_hat1;
   x0=b;

W = @(x) GNLM(x,x,patchRad,searchRad,hh);   % Linear NLM denoiser
% W= @(x) BM3D2(0, x, sigma, 'np', 0 );
% W= @(x) BM3DDEB2(x, sigma, h);

% W = @(x) guidedfilter(x,x,r0,eps);   % 
% [x_hat,converged,iters,errors,psnrs] = ...
%         pnpISTA2(x0,gradf,W,delta(i),x_orig,tol,maxiters);  % Run algorithm
[x_hat,converged,iters,errors,psnrs] = ...
        pnpISTA4(x0,gradf,W,delta(i),x_orig,tol,maxiters);  % Run algorithm
Dpsnrs(i,:)=psnrs;
lendsd=length(psnrs);
Dpsnrs(i,1:lendsd)=psnrs;
lendsd=length(errors);
Derrors(i,1:lendsd)=errors;

currentFile2 = sprintf('psnr_impaintingGUIDEDout2%d',delta(i))
folder_name2 = ['C:\Users\47058\Desktop\2020pnp\实验数据\' currentFile2  '.png']
imwrite(x_hat, folder_name2);

figure; imshow(x_hat); title('Output');
drawnow;
figure;
subplot(1,2,1); plot(errors,'LineWidth',2);
grid on; axis tight; xlabel('Iteration'); ylabel('Error');
title('Error');
subplot(1,2,2); plot(psnrs,'LineWidth',2);
grid on; axis tight; xlabel('Iteration'); ylabel('PSNR');
title('PSNR');


end

for i=1:1:length(delta)
    x0=b;
W = @(x) GNLM(x,x0,patchRad,searchRad,hh);   % Linear NLM denoiser
% W = @(x) guidedfilter(x0,x,r0,eps);   % 
% W= @(x) staticBM3DDEB2( x, sigma, h, 'np',x0 );
% W= @(x) staticBM3D( 0,x, sigma, 'np', 0 ,x0);
[x_hat,converged,iters,errors,psnrs] = ...
        pnpISTA4(x0,gradf,W,delta(i),x_orig,tol,maxiters);  % Run algorithm
% [x_hat,converged,iters,errors,psnrs] = ...
%         pnpISTA2(x0,gradf,W,delta(i),x_orig,tol,maxiters);  % Run algorithm
% Dpsnrs(i,:)=psnrs;
lendsd=length(psnrs);
Spsnrs(i,1:lendsd)=psnrs;

lendsd=length(errors);
Serrors(i,1:lendsd)=errors;
currentFile2 = sprintf('psnr_impaintingGUIDEDout0%d',delta(i))
folder_name2 = ['C:\Users\47058\Desktop\2020pnp\实验数据\' currentFile2  '.png']

imwrite(x_hat, folder_name2);


figure; imshow(x_hat); title('Output');
drawnow;
fig=figure;
subplot(1,2,1); plot(errors,'LineWidth',2);
grid on; axis tight; xlabel('Iteration'); ylabel('Error');
title('Error');
subplot(1,2,2); plot(psnrs,'LineWidth',2);
grid on; axis tight; xlabel('Iteration'); ylabel('PSNR');
title('PSNR');

drawnow;
end



%
function D=HtH(x)
h = fspecial('gaussian',[9 9],1);
H=h;
D=cconv2_by_fft2(x,H,0,[]);
D=cconv2_by_fft2(D,H,0,[]);
end
