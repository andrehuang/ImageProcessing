% File Name : main.m
function u = main(image_path , lambda_weight , ...
    kernel_size , gaussian_sigma , result_image_path, mu, tol)
if nargin <1 
    image_path="image/barbara.jpg" ; % 输入图片名
%     image_path = 'image/lena.jpg';
end 
if nargin <2 
    lambda_weight =2;  % lambda
end 
if nargin <3 
    kernel_size =15; % 卷积核大小
end 
if nargin <4 
    gaussian_sigma = 1 ; % 高斯核方差
end 
if nargin <5 
    result_image_path = "image/barbara_result.jpg"; % 输出图片名
%     result_image_path = 'image/lena_result.jpg';  
end 
if nargin < 6
    mu = 50;  % hyperparameter in ADMM in TVblur function；可对不同图像适当调整
end
if nargin < 7
    tol = 1e-4;  % stopping tolerance for ADMM iteration in TVblur function；可适当调整
end


kernel = fspecial('gaussian' , [kernel_size, kernel_size], gaussian_sigma );
img = im2double(imread(image_path));
% s = size(img); 
f = rgb2gray(img); % 读取图片
f = f/max(f(:));
figure;
imshow(f)

u = TVdeblur(f, kernel, lambda_weight, mu, tol); % 自己实现这个函数
h = figure;
imshow(u)
print(h, result_image_path , '-dpng') % 输出存储处理结果
end
