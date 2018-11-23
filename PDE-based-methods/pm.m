% File Name : pm.m
% can keep edges
function u = pm (image_path , T , K, result_image_path, fchoice)
if nargin <1 
%     image_path="image/barbara.jpg" ; % 输入图片名
    image_path = 'image/lena2.jpg';
end 
if nargin <2 
    T =0.3;  % Terminated time
end 
if nargin <3 
    K =1e-1;  
end 
if nargin <4
%     result_image_path = "image/barbara_pm_result.jpg"; % 输出图片名
    result_image_path = 'image/lena2_pm_result.jpg';
end 
if nargin <5
    fchoice = 2;
end

img = im2double(imread(image_path));
% [m,n] = size(img); 
f = rgb2gray(img); % 读取图片
% [m,n] = size(f);
f = f/max(f(:));
figure;
imshow(f)

% create noised picture
sigma1 = max(f(:))/20; % sd of noise
u0 = f + randn(size(f))*sigma1;
figure;
imshow(u0)

if fchoice == 1
    c = @(s) 1./(1+s./K);
end
if fchoice == 2
    c = @(s) 1./( 1+(s/K).^2 );
end

%     function amplitude = Amp(u)
%         s1= size(u);
%         amplitude =  (max(abs(grad(u)))).^2;
%         amplitude = reshape(amplitude, s1);
%     end
%     function delta = Delta(f)
%         [g1,g2] = gradient(f);
%         g1 = Amp(f) .* g1;
%         g2 = Amp(f) .* g2;
%         delta = divergence(g1, g2);
%     end

Amp = @(u) (abs(grad(u))).^2;
Delta = @(f) div(c(Amp(f)) .* grad(f));

tau = 0.01; % step size
niter = ceil(T/tau);

u = u0;
for i = 1:niter
    u = u + tau*Delta(u);
    psnr(u,f)
end
% u = TVdeblur(f, kernel, lambda_weight); % 自己实现这个函数
h = figure;
imshow(u)
print(h, result_image_path , '-dpng') % 输出存储处理结果
end
