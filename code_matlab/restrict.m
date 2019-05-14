function [uCoarse] = restrict(u)
% RESTRICT 

U = u;  

%Pad u with zeros to make the calculation simpler.
u = padarray(u,[1,1],0);

%Set uCoarse according to a linear combination of points from u.
uCoarse = .5*u(2:2:end-1,2:2:end-1) + .125*(u(1:2:end-2,2:2:end-1)...
                                          + u(3:2:end,2:2:end-1)...
                                          + u(2:2:end-1,1:2:end-2)...
                                          + u(2:2:end-1,3:2:end));

%Restrict u the old way by siply deleting every other point.
u = U(1:2:end,1:2:end);

%Reset the boundary values of uCoase with the correct boundary values.
uCoarse(:,1) = u(:,1);
uCoarse(:,end) = u(:,end);
uCoarse(1,:) = u(1,:);
uCoarse(end,:) = u(end,:);

end