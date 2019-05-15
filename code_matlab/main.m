close all
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
nValNum = length(minN:maxN);

stats = zeros(nValNum,maxN,2);
% hVec = 2.^(-(minN+1:maxN+1)+1);
hVec = 2.^(-(minN+1:maxN+1));

% errCell is for storing the error matrices at specific "frames"
% (iterations). frameCell is for remembering which "frames" those were.

errCell = cell(nValNum,maxN);
frameCell = cell(nValNum,maxN);

% exactSolCell stores g (the exact solution) evaluated on each mesh.

exactSolCell = cell(nValNum,1);
basesNum = 2;
iterVec = [5 50];

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
    exactSolCell{i-minN+1} = G;
    A = zeros(n);
    u0 = init(F,g,n,h,X,Y);

    u0(:,1) = g(X(:,1),Y(:,1));
    u0(:,n) = g(X(:,n),Y(:,n));
    u0(1,:) = g(X(1,:),Y(1,:));
    u0(n,:) = g(X(n,:),Y(n,:)); 

    N = n;

    for j = 1:4

        % The last argument is 0 because we don't (read: can't) use the mex
        % file, and that argument is what turns it on or off.
        [u,resMat,errMat,time,count] = looper2(F,g,n,N,j,2*iterVec,h,u0,xa,xb,ya,yb,tol,0);
        stats(i-minN+1,j,1) = norm(errMat(:,:,end),inf);
        stats(i-minN+1,j,2) = count;

        % Choose up to 4 almost evenly spaced iteration numbers and store them. 
        errFrames = unique(ceil(linspace(1,count,4)));
        frameCell{i-minN+1,j} = errFrames;
        
        % Store the error matrices at the chosen iterations.
        errCell{i-minN+1,j} = errMat(:,:,errFrames);
        
    end

end
    
% Save everything, just for safety.

save('stats.mat','stats');
save('errCell.mat','errCell');
save('exactSolCell.mat','exactSolCell');
   
%% Error and iteration number plots

legendStrs = {'One level','Two levels','Three levels','Four levels',...
    'Five levels','Six levels','Seven levels','Eight levels',...
    'Nine levels','Ten levels'};

% Error figure
errFig = figure;
plot(hVec,stats(:,1,1), 'o-')
xlabel('h')
ylabel('Error')
title('Error vs. h for all depths of recursion')
axis tight
% set(gca, 'XDir','reverse') %reverses the order of h
saveas(errFig,'errFig.fig')

% Iteration figure
countFig = figure;
semilogy(hVec,stats(:,:,2),'o-');
legend(legendStrs(1:4));
xlabel('h')
ylabel('Iterations')
title(sprintf('Number of iterations vs. h for %d depths of recursion',4))
axis tight
% set(gca, 'XDir','reverse')
saveas(countFig,'countFig.fig')

%% Error surface plots

errorDir = 'error_surfs';
mkdir(errorDir);
exactSolDir = 'exact_sol_surfs';
mkdir(exactSolDir);

for i = 1:nValNum
    
    nValDir = sprintf('%s/N_%d',errorDir,1/hVec(i)+1);
    mkdir(nValDir)
   
    for j = 1:maxN
        if ~isempty(errCell{i,j}) 
            depthDir = sprintf('%s/depth_%d',nValDir,j);
            mkdir(depthDir)
            subplotNum = size(errCell{i,j},3);
            err = errCell{i,j};
            
            for k = 1:subplotNum
                fig1 = figure;
                surf(abs(err(:,:,k)),'linestyle','none');
                title(sprintf('h = %f, depth = %d levels, iteration = %d',...
                    hVec(i),j,frameCell{i,j}(k)));
                zlim([0 norm(err(:),inf)]);
                saveas(fig1,sprintf('%s/count_%d.fig',depthDir,frameCell{i,j}(k)));
            end
        end
    end
    
    fig2 = figure;
    surf(exactSolCell{i},'linestyle','none');
    title(sprintf('Exact solution evaluated for h = %f',hVec(i)));
    saveas(fig2,sprintf('%s/N_%d.fig',exactSolDir,1/hVec(i)+1));
end

%% Error Surface Plots 

% error at N = 65, depth 2
err = errCell{2,2}

fig1 = figure;
surf(linspace(0,1,65),linspace(0,1,65),err(:,:,1));
title('error evaluated at N=65, depth_2, count 1')
saveas(fig1,sprintf('error_65_d2_c1.fig'))

fig2 = figure;
surf(linspace(0,1,65),linspace(0,1,65),err(:,:,2));
title('error evaluated at N=65, depth_2, count 2')
saveas(fig2,sprintf('error_65_d2_c2.fig'))

fig3 = figure;
surf(linspace(0,1,65),linspace(0,1,65),err(:,:,3));
title('error evaluated at N=65, depth_2, count 3')
saveas(fig3,sprintf('error_65_d2_c3.fig'))

fig4 = figure;
surf(linspace(0,1,65),linspace(0,1,65),err(:,:,4));
title('error evaluated at N=65, depth_2, count 4')
saveas(fig4,sprintf('error_65_d2_c4.fig'))