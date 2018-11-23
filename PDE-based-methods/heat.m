% File Name : heat.m
function u = heat(image_path , T , result_image_path)
if nargin <1 
    image_path="image/barbara.jpg" ; % 输入图片名
%     image_path = 'image/lena2.jpg';
end 
if nargin <2 
    T =5e-6;  % Terminated time
end 
if nargin <3
    result_image_path = "image/barbara_heat_result.jpg"; % 输出图片名
%     result_image_path = 'image/lena2_heat_result.jpg';
end 

img = im2double(imread(image_path));
s = size(img); 
n = s(1);
f = rgb2gray(img); % 读取图片
f = f/max(f(:));
figure;
imshow(f)

% create noised picture
sigma1 = max(f(:))/10; % sd of noise
u0 = f + randn(size(f))*sigma1;
psnr(u0,f)
figure;
imshow(u0)


% 
h = 1/n;
nabla = fspecial('laplacian',0);
Delta = @(f) 1/h^2 * imfilter(f, nabla, 'circular');
tau = 1e-7;
niter = ceil(T/tau);
u = u0;
for i = 1:niter
    u = u + tau*Delta(u);
    psnr(u,f)
end
% psnr(u,f)

% u = TVdeblur(f, kernel, lambda_weight); % 自己实现这个函数
h = figure;
imshow(u)
print(h, result_image_path , '-dpng') % 输出存储处理结果
end