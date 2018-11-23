function g = grad(U) % forward scheme
%         [FX, FY] = gradient(x);
        [mx, nx] = size(U);
        g = zeros([2,mx,nx]);
        
        FX = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
        FY = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
            
        g(1,:,:) = FX;
        g(2,:,:) = FY;
end