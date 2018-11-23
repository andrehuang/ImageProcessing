% File Name : shockfilters.m
function u = shockfilters(image_path , T , result_image_path, Lchoice)
if nargin <1 
%     image_path="image/fig1.png" ; % 输入图片名
    image_path = 'image/lena2.jpg';
%     image_path = 'image/yosemite.png';
end 
if nargin <2 
    T =50;  % Terminated time
end 
if nargin <3
%     result_image_path = "image/fig1_shock_result.png"; % 输出图片名
    result_image_path = 'image/lena2_shock_result.jpg';
%     result_image_path = 'image/yosemite_shock_result.png';

end 
if nargin <4
    Lchoice=2;
end


noise = true;
tau = 0.1; % step size

img = im2double(imread(image_path));
f = rgb2gray(img); 
s = size(f); 
n = s(1);
f = f/max(f(:));
figure;
imshow(f)

% blur kernel
gaussian_sigma = 2;
kernel_size = 10;
kernel = fspecial('gaussian' , [kernel_size, kernel_size], gaussian_sigma );

if noise
    sigma1 = max(f(:))/200; % sd of noise
    u0 = imfilter(f,kernel, 'circular') + randn(size(f))*sigma1;
else
    u0 = imfilter(f,kernel, 'circular');
end

figure;
imshow(u0)
% 
% h = 1/n;
% GradAbs = @(f) reshape(sqrt(sum(grad(f).^2, 1)),s);

    function [Lf, a_grad_I] = L1(I)
        nabla = fspecial('laplacian',0);
        Lf = imfilter(f, nabla, 'circular');
        % center difference grad
        [nx, ny] = size(I);
        I_mx = I-I(:,[1 1:ny-1]);
        I_px = I(:,[2:ny ny])-I;
   
        I_my = I-I([1 1:nx-1],:);
        I_py = I([2:nx nx],:)-I;
        % minmod operator
        Dx = min(abs(I_mx),abs(I_px));
        ind=find(I_mx.*I_px < 0); Dx(ind)=zeros(size(ind));
        Dy = min(abs(I_my),abs(I_py));
        ind=find(I_my.*I_py < 0); Dy(ind)=zeros(size(ind));
        
        a_grad_I = sqrt(Dx.^2+Dy.^2);
    end

    function [I_nn, a_grad_I] = L2(I)
        % center difference grad
        [nx, ny] = size(I);
        I_mx = I-I(:,[1 1:ny-1]);
        I_px = I(:,[2:ny ny])-I;
   
        I_my = I-I([1 1:nx-1],:);
        I_py = I([2:nx nx],:)-I;
   
        I_x = (I_mx+I_px)/2;
        I_y = (I_my+I_py)/2;
        % minmod operator
        Dx = min(abs(I_mx),abs(I_px));
        ind=find(I_mx.*I_px < 0); Dx(ind)=zeros(size(ind));
        Dy = min(abs(I_my),abs(I_py));
        ind=find(I_my.*I_py < 0); Dy(ind)=zeros(size(ind));
        
        a_grad_I = sqrt(Dx.^2+Dy.^2);
        % center difference second-order direvative
        I_xx = I(:,[2:ny ny])+I(:,[1 1:ny-1])-2*I;
        I_yy = I([2:nx nx],:)+I([1 1:nx-1],:)-2*I;
        I_xy = (I_x([2:nx nx],:)-I_x([1 1:nx-1],:))/2;
        
        

%         FXX = [diff(FX,1,2), FX(:,1,:) - FX(:,end,:)];
%         FYY = [diff(FY,1,1); FY(1,:,:) - FY(end,:,:)];
%         FXY = [diff(FX,1,1); FX(1,:,:) - FX(end,:,:)];
%         Lf  = 1./(GradAbs(U)).^2 .*(FX.*FX .* FXX + 2.*FX.*FY.*FXY + FY.*FY .* FYY);
        I_nn = 1./(a_grad_I).^2 .* (I_x.*I_x.*I_xx + 2*I_x.*I_y.*I_xy + I_y.*I_y.*I_yy);
%         dl = 1e-8;
%         I_nn = I_xx.*abs(I_x).^2 + 2*I_xy.*I_x.*I_y + I_yy.*abs(I_y).^2;
%         I_nn = I_nn./(abs(I_x).^2+abs(I_y).^2+dl);
%         a2_grad_I = (abs(I_x)+abs(I_y)); % second order abs grad
%         ind = find(a2_grad_I==0);  % zero gradient n,e not defined
%         I_nn(ind)=I_xx(ind); 
    end

if Lchoice == 1
    L = @L1;
end
if Lchoice == 2
    L = @L2;
end

% Delta = @(f) -GradAbs(f) .* sign(L(f));

niter = ceil(T/tau);
u = u0;
for i = 1:niter
    [Lf, a_grad_I] = L(u);
    delta = -a_grad_I .* sign(Lf);
    u = u + tau*delta;
    psnr(u,f);
    ssim(u,f);
end


% u = TVdeblur(f, kernel, lambda_weight); % 自己实现这个函数
h = figure;
imshow(u)
print(h, result_image_path , '-dpng') % 输出存储处理结果
end
