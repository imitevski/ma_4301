function [u,resMat,err] = FAS_V2(F,g,n,N,levels,iterVec,h,u0,xa,xb,ya,yb,count,mex)
%FAS_V2.m implements the Full Approximation Scheme. 

direction = 'down';
plotFigs = 0;

%Do the initial Gauss-Seidel iteration.
[u,resMat] = GaussSeidel(F,g,iterVec,h,u0,xa,xb,ya,yb,0,mex);
res = norm(resMat(:),inf);

if plotFigs == 1
    fig1 = figure(18);
    surf(resMat,'linestyle','none')
    title(sprintf('n = %d, res = %f, count = %d',n,res,count))
    drawnow
end

%Go down one level of granularity, setting a new n and h.
n = (n+1)/2;
h = 2*h;

%Set the coarse versions of u, F, and the residual matrix.
v = restrict(u);
resCoarse = restrict(resMat);

%Evaluate the MA equation on the coarse grid, pad it with zeros so that
%it's the same size as resCoarse, and add them together.
A = padarray(A_solver(v,h,2),[1,1],0) + resCoarse;

n

%If we're on the lowest level, run the Gauss-Seidel calculation and set
%vNew as the output. Otherwise, set vNew as the recursive output of FAS_V.
if floor(log2(N)) - floor(log2(n)) ~= levels

    vNew = FAS_V2(A,g,n,N,levels,iterVec,h,v,xa,xb,ya,yb,count,mex);
    
else
    
    vNew = GaussSeidel(A,g,iterVec,h,v,xa,xb,ya,yb,1,mex);
    
end

direction = 'up';

%Calculate the coarse error matrix.
eCoarse = vNew - v;

%Interpolate the coarse error matrix to the fine level.
eFine = enhance(eCoarse);

%Add the error matrix to the original matrix u.
u(2:end-1,2:end-1) = u(2:end-1,2:end-1) + eFine(2:end-1,2:end-1);

%Come back up to the fine level.
n = 2*n-1;
h = h/2;

%Do a few more Gauss-Seidel iterations on u and calculate the residual 
%matrix.
[u,resMat] = GaussSeidel(F,g,iterVec,h,u,xa,xb,ya,yb,0,mex);
res = norm(resMat(:),inf);

n

if plotFigs == 1
    fig2 = figure(36);
    surf(resMat,'linestyle','none')
    title(sprintf('n = %d, res = %f, count = %d',n,res,count))
    drawnow
end

%If we're on the finest level, calculate the error. Otherwise, set it to
%zero so as to avoid a 'not enough output arguments' error.
if n == N
    [X,Y] = meshgrid(xa:h:xb,ya:h:yb);
    G = g(X,Y);
    err = G-u;
else
    err = 0;
end

end