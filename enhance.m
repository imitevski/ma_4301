function [uFine] = enhance(u)
% ENHANCE does interpolation when moving from coarse to fine level

n = 2*length(u)-1;
uFine = zeros(n);

uFine(1:2:end,1:2:end) = u(:,:);

uFine(2:2:end-1,1:2:end) = .5*(u(1:end-1,:) + u(2:end,:));
uFine(1:2:end,2:2:end-1) = .5*(u(:,1:end-1) + u(:,2:end));
uFine(2:2:end-1,2:2:end-1) = .25*(u(1:end-1,1:end-1) + u(2:end,1:end-1)...
                                + u(1:end-1,2:end) + u(2:end,2:end));
                            
end