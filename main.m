% close all
clear all

% The solver is executed by running this file

% f1 and f2 are different depending on which basis we use
f1 = @(x,y) (1 + x.^2).*(1 + y.^2).*exp(x.^2 + y.^2);
f2 = @(x,y) exp(x.^2 + y.^2).*(1+.5*(x+y).^2).*(1+.5*(y-x).^2);
% f = @(x,y) (1 + x.^2 + y.^2).*exp(x.^2 + y.^2);
g = @(x,y) exp((x.^2 + y.^2)/2);

%Set the minimum and maximum values of n for the solution matrix.
minN = 4;
maxN = 7;
stats = zeros(length(minN:maxN),maxN+1,3);
stats(:,1,:) = repmat((2.^(minN+1:maxN+1)+1)',[1 1 3]);

basesNum = 2;
iterVec = [5 50];
error_1 = zeros(minN:maxN)

% k == 0 if we don't want to use the mex file, k == 1 if we do.
for k = 0:0
    
    for i = minN:maxN

        i

        n = 2^(i+1) + 1;
        h = 1/(n-1);
        xa = 0; xb = 1; ya = 0; yb = 1; tol = h^2/10;
        [X,Y] = meshgrid(xa:h:xb,ya:h:yb);

        if basesNum == 1

            F = f1(X,Y);

        else

            F1 = f1(X,Y);
            F2 = f2(X,Y);

            F = min(F1,F2);

        end

        G = g(X,Y);

        A = zeros(n);

        u0 = init(F,g,n,h,X,Y);

        u0(:,1) = g(X(:,1),Y(:,1));
        u0(:,n) = g(X(:,n),Y(:,n));
        u0(1,:) = g(X(1,:),Y(1,:));
        u0(n,:) = g(X(n,:),Y(n,:)); 

    %     %Calculate u_xx*u_yy, and leave only the nonnegative elements.
    %     u_xy = max(1./h.^2.*(u0(1:end-2,2:end-1)+u0(3:end,2:end-1)-2.*u0(2:end-1,2:end-1)),0)...
    %        .*max(1./h.^2.*(u0(2:end-1,1:end-2)+u0(2:end-1,3:end)-2.*u0(2:end-1,2:end-1)),0);
    %    
    %     %Calculate u_vv*u_v(perp)v(perp), and leave only the nonnegative elements.
    %     u_vv = max(1./(2.*h.^2).*(u0(3:end,1:end-2)+u0(1:end-2,3:end)-2.*u0(2:end-1,2:end-1)),0)...
    %            .*max(1./(2.*h.^2).*(u0(3:end,3:end)+u0(1:end-2,1:end-2)-2.*u0(2:end-1,2:end-1)),0);
    %        
    %     mask = u_xy <= u_vv;
    %     
    %     pcolor(mask)

        N = n;

        for j = 1:i

            [u,resMat,error,time,count] = looper2(F,g,n,N,j,2*iterVec,h,u0,xa,xb,ya,yb,tol,k); 
%             stats(i,j+1,k+1) = count;
            stats(i,j+1,k+1) = norm(error,inf);
%             print("The error is ", error)
            save('stats.mat','stats');
%         error_1(i) = norm(error,inf)
        
        end

    end
    
end
save('stats.mat','stats');
   
fig = figure;
hold on

for i = 1:maxN
    
    semilogy(stats(i:end,1),stats(i:end,i+1))
    
end

legend('One level','Two levels','Three levels','Four levels',...
    'Five levels','Six levels','Seven levels','Eight levels');
xlabel('N')
ylabel('log(Time)')
title('Time performance vs. N for eight depths of recursion')
    
saveas(fig,'time_stats.jpg')
