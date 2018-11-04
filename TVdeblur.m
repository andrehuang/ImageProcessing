function u = TVdeblur(u_original, kernel, lambda_weight, mu, tol)
% generate blurred and noised picture f
sigma1 = max(max(u_original))/100;
f = imfilter(u_original, kernel) + randn(size(u_original))*sigma1;
[m,n] = size(u_original);
figure;
imshow(f)
psnr1 = psnr(f,u_original)
ssim1 = ssim(f,u_original)
% =======   use admm to solve the TV deblur problem  =======
% hyper-parameters
% mu = 50;
% tol = 1e-4; % depend on the picture, 1e-4, 1e-5, 1e-6; it should be bigger when image size is larger
delta = 1;
% variable initialization
d = rand([2,m,n]);
b = rand([2,m,n]);
u = f;

% =============  function and operators      =========
% operator = @(x) reshape(imfilter(imfilter(x, kernel), kernel') + mu*del2(x),[m*n,1]);
    function DxyT = div(x)  
        % Transpose of the forward finite difference operator: 
        %    the divergence of the forward finite difference operator
        % input is a 2*m*n array, the gradient vector of an image
        [~,m1,n1] = size(x);
        X = reshape(x(1,:,:),[m1,n1]);
        Y = reshape(x(2,:,:),[m1,n1]);
        
        [rows, cols] = size(X);
        dxt = zeros(rows, cols);
        dxt(:,1:cols-1) = X(:,1:cols-1)-X(:,2:cols);
        dxt(:,cols) = X(:,cols)-X(:,1);
        

        dyt = zeros(rows, cols);
        dyt(1:rows-1,:) = Y(1:rows-1,:)-Y(2:rows,:);
        dyt(rows,:) = Y(rows,:)-Y(1,:);

        
%         dxt = [X(:,end) - X(:, 1), -diff(X,1,2)];
%         dyt = [Y(end,:) - Y(1, :); -diff(Y,1,1)];

        DxyT = dxt + dyt;
    end

    function g = grad(U) % forward scheme
%         [FX, FY] = gradient(x);
        [mx, nx] = size(U);
        g = zeros([2,mx,nx]);
        
        FX = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
        FY = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
%         FX = [diff(U,1,2), zeros(
            
        g(1,:,:) = FX;
        g(2,:,:) = FY;
    end

    function a = three_norm(arr)
        [~,y,z] = size(arr);
%         a = zeros([y,z]);
        a = arr(1,:,:) + arr(1,:,:);
        a = reshape(a, [y,z]);
        a = norm(a);
    end
    
    function z = shrink(x,r)
        z = sign(x).*max(abs(x)-r,0);
    end

A   = psf2otf(kernel,[m,n]); %computes the FFT of the point-spread kernel
laplaceF   = fspecial('laplacian',0);  
laplace = psf2otf(laplaceF, [m n]);
lhs     = mu*conj(A).*A - lambda_weight*laplace;
Atf     = imfilter(f, kernel', 'circular');

% ============  begin admm iteration  ==============
while three_norm(grad(u)-d)/norm(f) > tol
    rhs = mu*Atf + lambda_weight*div(d-b); 
    u = fft2(rhs)./lhs;
    u = real(ifft2(u)); 
%     norm(u) 
    d = shrink(grad(u) + b, lambda_weight/mu);
    b = b + delta*(grad(u)-d);
%     three_norm(grad(u)-d)/norm(f);

end
psnr2 = psnr(u,u_original)
ssim2 = ssim(u,u_original)


end