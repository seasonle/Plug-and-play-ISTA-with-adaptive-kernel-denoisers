clear; close all; clc;
addpath('./utils');

%% Forward model
% x_orig = im2double(imread('./house.tif')); % Ground-truth
% x_orig = im2double(imread('3_peppers256.png')); % Ground-trut
% x_orig = x_orig(2:256,2:256);
% x_orig = im2double(imread('6_boat.png')); % Ground-trut
x_orig = im2double(imread('2_house.png')); % Ground-trut

% x_orig =x_orig(1:100,1:100) ; % Ground-trut
% x_orig = im2double(imread('house.tif')); % Ground-trut
%  x_orig = im2double(imread('4_Lena512.png')); % Ground-trut
% x_orig = im2double(imread('house512.tif')); % Ground-trut
% x_orig = im2double(imread('boats.tif')); % Ground-trut
% x_orig = im2double(imread('7_hill.png')); % Ground-trut
% x_orig = im2double(imread('4_Lena256.tif')); % Ground-trut
% x_orig = im2double(imread('1_Cameraman256.png')); % Ground-trut

%  x_orig = im2double(imread('4_Lena512.png')); % Ground-trut
% x_orig = im2double(imread('barbara.tif')); % Ground-trut
% x_orig = im2double(imread('peppers.tif')); % Ground-trut
% r = 0.7;                % Fraction of missing pixels0.7
 r = 0.8;   
sigma_n = 10/255;       % Noise standard deviation (on a scale of [0,1])

[rr,cc] = size(x_orig);
P = double(rand(rr,cc)>=r);     % Decimation operator
b = P.*x_orig;
b(P==1) = b(P==1) + sigma_n*randn(nnz(P),1);    % Observed image
b(b>1) = 1; b(b<0) = 0;
gradf = @(x) P.*(x-b);      % Gradient of data fidelity term (f)

imshow(b); title('Input (decimated) image');
pause(0.1); drawnow;

%% Algorithm parameters
searchRad = 10;             % Search window radius in NLM
patchRad = 3;               % Patch radius in NLM
% sigma=10;%10

tol = 1e-5;                 % Tolerance for termination of algorithm
maxiters = 300;             % Max. no. of iterations
% delta = [0.3,0.5,0.7,0.9,1.0,1.2,1.5,1.7,1.9,2.0,2.1];  
% delta = [0.3,0.5,0.7,0.9,1.0,1.5,1.9,2.0,2.1];  
delta = [0.6,0.9,1.2,1.6,1.8];
% delta = 1.8;                % Step-size in PnP-ISTA
% delta = [1.6];                % Step-size in PnP-ISTA
% delta = [0.3,0.6,0.9,1.0,1.2,1.5,1.7,1.9,2.0,2.1];
% delta = [0.9];
%% Main algorithm
x0 = initInpainting(b,P==0,5);  % Initial point to start the iterations
% W= @(x) staticBM3D(0, x, sigma, 'np', 0 ,x0);
figure; imshow(x0); title('Initialization');
Pxb=psnr(b,x_orig)
Px0=psnr(x0,x_orig)
%
h=12/255;%0.7时为15
h1=12/255;
  h2=12/255;
% W1 = @(x) JNLM(x,x0,patchRad,searchRad,h1);
W1 = @(x) GNLM(x,x0,patchRad,searchRad,h1);

Sx0NLM=cal_ssim(x0.*255,x_orig.*255,0,0)
Xnlm=W1(x0);
Pxnlm=psnr(Xnlm,x_orig)
SfNLM=cal_ssim(Xnlm.*255,x_orig.*255,0,0)


% currentFile0 = sprintf('psnr_impaintingBM3D0input%d');
% folder_name0 = ['C:\Users\xingle\Desktop\2020pnp\实验数据\BM3D\BM3D0\' currentFile0  '.png'];
% currentFile00 = sprintf('psnr_impaintingBM3D0x0%d');
% folder_name00 = ['C:\Users\xingle\Desktop\2020pnp\实验数据\BM3D\BM3D0\' currentFile00  '.png'];
% currentFile01 = sprintf('psnr_impaintingoriginalnlminput%d');
% folder_name01 = ['C:\Users\xingle\Desktop\2020pnp\实验数据\BM3D\BM3D0\' currentFile01  '.png'];


currentFile0 = sprintf('psnr_impaintingBM3D0input%d');
folder_name0 = ['C:\Users\47058\Desktop\2020pnp\实验数据\BM3D\' currentFile0  '.png'];
currentFile00 = sprintf('psnr_impaintingBM3D0x0%d');
folder_name00 = ['C:\Users\47058\Desktop\2020pnp\实验数据\BM3D\' currentFile00  '.png'];
currentFile01 = sprintf('psnr_impaintingoriginalnlminput%d');
folder_name01 = ['C:\Users\47058\Desktop\2020pnp\实验数据\BM3D\' currentFile01  '.png'];



imwrite(b, folder_name0);
imwrite(x0, folder_name00);
imwrite(Xnlm, folder_name01);
% DL=length(delta);
% Dpsnrs=zeros(DL,maxiters);
Spsnrs=zeros(length(delta),maxiters);
Dpsnrs=zeros(length(delta),maxiters);


Sfinalssim=zeros(length(delta),1);
Dfinalssim=zeros(length(delta),1);
% Sfinalpsnr=zeros(length(delta),1);
% Dfinalpsnr=zeros(length(delta),1);




  
%   maxiters=100;
W1 = @(x) GNLM(x,x0,patchRad,searchRad,h2); 
% W1 = @(x) JNLM(x,x0,patchRad,searchRad,h);  

for i=1:1:length(delta)
  
[x_hat1,converged,iters,errors,psnrs1] = ...
        pnpISTA4(x0,gradf,W1,delta(i),x_orig,tol,10);  % Run algorithm
    
    
x00=x_hat1;
% W = @(x) JNLM(x,x,patchRad,searchRad,h);   % Linear NLM denoiser

W = @(x) GNLM(x,x,patchRad,searchRad,h2); 
% W= @(x) BM3D2(0, x, sigma, 'np', 0 );

% W = @(x) GNLM(x,x,patchRad,searchRad,h); 
% W = @(x) guidedfilter(x,x,r0,eps);
[x_hat,converged,iters,errors,psnrs] = ...
        pnpISTA4(x00,gradf,W,delta(i),x_orig,tol,maxiters);  % Run algorithm
% Dpsnrs(i,:)=psnrs;
lendsd=length(psnrs);
Dpsnrs(i,1:10)=psnrs1;
Dpsnrs(i,11:lendsd)=psnrs(1:lendsd-10);

% lendsd=length(psnrs);
% Dpsnrs(i,1:lendsd)=psnrs;
Dfinalssim(i)=cal_ssim(x_hat.*255,x_orig.*255,0,0)
fprintf('\nNo. of iterations completed = %d\n',iters);
currentFile1 = sprintf('psnr_impaintingBM3D2%d',delta(i))
folder_name1 = ['C:\Users\47058\Desktop\2020pnp\实验数据\BM3D2\' currentFile1  '.mat']
save(folder_name1,'psnrs')%save 

% currentFile2 = sprintf('psnr_impaintingnlm%d',delta(i))
currentFile2 = sprintf('psnr_impaintingBM3D2out%d',delta(i))
folder_name2 = ['C:\Users\47058\Desktop\2020pnp\实验数据\BM3D2\' currentFile2  '.png']

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

currentFile3 = sprintf('psnr_impaintingBM3D2%d',delta(i))

folder_name3 = ['C:\Users\47058\Desktop\2020pnp\实验数据\BM3D2\' currentFile3  '.png']
saveas(gcf,folder_name3,'png')
drawnow;
end






% 
for i=1:1:length(delta)
% W = @(x) GNLM(x,x0,patchRad,searchRad,h); 
% W = @(x) guidedfilter(x0,x,r0,eps);
% W= @(x) staticBM3D(0, x, sigma, 'np', 0 ,x0);
% W2 = @(x) JNLM(x,x0,patchRad,searchRad,h);   % Linear NLM denoiser
 W2 = @(x) GNLM(x,x0,patchRad,searchRad,h); 
[x_hat,converged,iters,errors,psnrs] = ...
        pnpISTA4(x0,gradf,W2,delta(i),x_orig,tol,maxiters);  % Run algorithm
% Dpsnrs(i,:)=psnrs;
lendsd=length(psnrs);
Spsnrs(i,1:lendsd)=psnrs;
Sfinalssim(i)=cal_ssim(x_hat.*255,x_orig.*255,0,0)
fprintf('\nNo. of iterations completed = %d\n',iters);
currentFile1 = sprintf('psnr_impaintingBM3D0%d',delta(i))
folder_name1 = ['C:\Users\47058\Desktop\2020pnp\实验数据\BM3D\' currentFile1  '.mat']
save(folder_name1,'psnrs')%save 

% currentFile2 = sprintf('psnr_impaintingnlm%d',delta(i))
currentFile2 = sprintf('psnr_impaintingBM3Dout0%d',delta(i))
folder_name2 = ['C:\Users\47058\Desktop\2020pnp\实验数据\BM3D\' currentFile2  '.png']

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
frame=getframe(fig);
imgplot=frame2im(frame);
currentFile3 = sprintf('psnr_impaintingBM3D0%d',delta(i))

folder_name3 = ['C:\Users\47058\Desktop\2020pnp\实验数据\BM3D\' currentFile3  '.png']
% imwrite(imgplot, folder_name3);
saveas(gcf,folder_name3,'png')
drawnow;
end




