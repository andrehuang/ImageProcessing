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