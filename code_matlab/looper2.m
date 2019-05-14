function [u,resMat,err,time,count] = looper2(F,g,n,N,levels,iterVec,h,u0,xa,xb,ya,yb,tol,mex)

%LOOPER2 is used to put stopping critera on the FAS_V2 iterations

tic

count = 0;

u = u0;
res = 1;
resLast = 1e2;

while res > tol && count < 15
    
    count = count + 1;
    
    [u,resMat,err] = FAS_V2(F,g,n,N,levels,iterVec,h,u,xa,xb,ya,yb,count,mex);
    res = norm(resMat(:),inf);
%     error = norm(err,inf)
    
    resLast = res;
        
end

time = toc;


end