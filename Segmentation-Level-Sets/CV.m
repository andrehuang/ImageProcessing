function PhiOut = CV(image_path, result_image_path)
if nargin <1 
    image_path="image/cortex.png" ; % 输入图片名
%     image_path = 'image/lena2.jpg';
end 
if nargin <2
    result_image_path = "image/kitty_GAC_result.jpg"; % 输出图片名
%     result_image_path = 'image/lena2_heat_result.jpg';
end 
img = im2double(imread(image_path));
s = size(img) ;
m = s(1);
n = s(2);
if length(s) >= 3
    f = rgb2gray(img); % 读取图片
else
    f = img;
end
I = f/max(f(:));

% [Y,X] = meshgrid(1:n,1:m);
% r = max(m,n)/3;
% c = [m n]/2;
% PhiIn = max( abs(X-c(1)), abs(Y-c(2)) ) - r;


PhiIn = ones(m,n)*(-1);
PhiIn(3:m-3, 3:n-3) = 1;
PhiIn(3,3:n-3) = 0;
PhiIn(m-3,3:n-3) = 0;
PhiIn(3:m-3,n-3) = 0;
PhiIn(3:m-3,3) = 0;

param.tau = 0.01;		%Time step.
param.kernel_size = 5;		
param.gaussian_sigma=2;
param.p = 1;
param.Iter = 300;		%Number of iterations
param.Rein_Iter = 10;
param.alpha = 1;
param.eps = 1e-3;
options.order = 1;

% PhiOut = Reinitial2D(PhiIn, 10);

kernel = fspecial('gaussian', param.kernel_size, param.gaussian_sigma);
smoothed_I = imfilter(I, kernel, 'replicate');
% figure;
% imshow(smoothed_I)
% [gradIx, gradIy] = gradient(smoothed_I);
% gradI =  sqrt(gradIx.^2 + gradIy.^2);
gradientI = grad(smoothed_I,options);
gradI = sqrt(sum(gradientI.^2,3));
% p = param.p;
% W = 1./ (param.eps+gradI.^p);
W = 1./(param.eps+gradI);
W = rescale(W,0,1);
% W = rescale(-gradI, 0.1, 1);
figure;
imshow(W);

% g = grad(smoothed_I,options);
% d0 = sqrt(sum(g.^2,3));
% W = 1./(param.eps+d0);
% [gradgx, gradsgy] = gradient(g);
gW = grad(W,options);
% laplaceF   = fspecial('laplacian',0);
Phi = PhiIn;
lambda=0.8;
c1 = 0.7;c2=0;
for i = 1: param.Iter
%     break
    gPhi = grad(Phi,options);
    absgPhi = max(param.eps, sqrt(sum(gPhi.^2,3)) );
    g = gPhi ./ repmat( absgPhi, [1 1 2] );
%     G =  W .* absgPhi .* div(g,options ) + ...
%         param.alpha * W .* absgPhi + sum(gW.*gPhi,3);
    G = absgPhi .* div( g,options ) - lambda*(I-c1).^2 + lambda*(I-c2).^2;

%     [gradPhix, gradPhiy] = gradient(Phi);
%     gradPhi = max(param.eps, sqrt(gradPhix.^2 + gradPhiy.^2));
% %     G1 = g .* gradPhi .* div(grad(Phi)./reshape(repmat(gradPhi,1,1,2), [2,m,n]));
%     G1 = g.*gradPhi.*imfilter(Phi./gradPhi, laplaceF,'replicate');
%     G2 = param.alpha * (max(g,0).*Delta_pos(Phi) + min(g,0).*Delta_neg(Phi));
%     G3 = max(gradgx,0).*deltax_neg(Phi) + min(gradgx,0).*deltax_pos(Phi) + max(gradgy,0).*deltay_neg(Phi) + min(gradgy,0).*deltay_pos(Phi); 
%     delta = G1 + G2 + G3;
    Phi = Phi + param.tau * G;
%     imshow(Phi)
%     contour(Phi, [0,0],'r-')
    imshow(I),hold on, contour(Phi, [0,0],'r-', 'LineWidth', 1.2), title([int2str(i),'/',int2str(param.Iter)]), drawnow,
	hold off

%     hold on;
%     imshow(I);
%     [~,h] = contour(Phi, [t t], 'r');
%     set(h, 'LineWidth', 2);
%     axis ij;
%     hold off;
    if mod(i, param.Rein_Iter) == 0
        Phi = Reinitial2D(Phi, 10);
    end
end

PhiOut = Phi;
    function d = deltax_pos(u)
        d = [diff(u,1,2), u(:,1,:) - u(:,end,:)];
    end

    function d = deltax_neg(u)
        [~,col] = size(u);
        d = [u(:,1:col-1)-u(:,2:col), u(:,col)-u(:,1)];
    end

    function d =deltay_pos(u)
        [row, ~] =size(u);
        d = [u(2:row,:)-u(1:row-1,:); u(1,:)-u(row,:)];
    end

    function d = deltay_neg(u)
        [row,~] = size(u);
        d = [u(1:row-1,:)-u(2:row,:);u(row,:)-u(1,:)];
    end

    function d = Delta_pos(u)
        d=max(deltax_neg(u),0).^2+min(deltax_pos(u),0).^2 + max(deltay_neg(u),0).^2 + min(deltay_pos(u),0).^2;
        d = sqrt(d);
    end

    function d = Delta_neg(u)
        d=min(deltax_neg(u),0).^2+max(deltax_pos(u),0).^2 + min(deltay_neg(u),0).^2 + max(deltay_pos(u),0).^2;
        d = sqrt(d);
    end

end
