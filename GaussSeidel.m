function [u,resMat,err] = GaussSeidel(F,g,iterVec,h,u0,xa,xb,ya,yb,coarse,mex)
% GaussSeidel does Gauss Seidel iterations

%Lambda is a constant used for calculating a weighted combination of the
%existing u and the updated u.
lambda = 0.05;

%Set the grid on which the exact solution g(x,y) will be applied.
[X,Y] = meshgrid(xa:h:xb,ya:h:yb);
G = g(X,Y);

%Initialize the matrix u.
u = u0;

%Set the right-hand side of the MA equation and the number of iterations 
%in the loop depending on whether the input matrix u is on the coarse or 
%the fine level.
if coarse == 1
    maxCount = iterVec(1);
else
    maxCount = iterVec(2);
end

for cuenta = 1:maxCount
    
    %mex is always 0 for these runs 
    if mex == 0
        uNew = notJacobi(u,F,h,2);
    else
        [m,n] = size(u);
        [uNew,hInv] = notNotGaussSeidel(u,F,h,m,n);
        uNew = uNew(2:end-1,2:end-1);
    end
    
%     figure(18)
%     u_Mat = notJacobi(u,F,h,2);
%     
%     [m,n] = size(u);
%     u_C = notNotGaussSeidel(u,F,h,m,n);
%     u_C = u_C(2:end-1,2:end-1);
    
%     surf(abs(u_Mat - u_C))
%     surf(u_C)
    
%     uNew = notJacobi(u,F,h,2);
    
    %Update u with a weighted sum of points from the already-extant u and
    %the newly calculated uNew.
    u(2:end-1,2:end-1) = lambda*u(2:end-1,2:end-1) + (1-lambda)*uNew;
    
    resMat = padarray(F(2:end-1,2:end-1) - A_solver(u,h,2),[1,1],0);
    
%     figure(18)
%     surf(resMat)
%     drawnow
    
end

%Subtract u from the exact solution to find the error matrix.
err = G - u;

%Pad the residual matrix with zeros to make it the right size again.
resMat = padarray(F(2:end-1,2:end-1) - A_solver(u,h,2),[1,1],0);

end