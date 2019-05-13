function [U] = A_solver(u,h,basesNum)
% A_SOLVER calculates the LHS of the equation using 2 basis and brings back
% the non-negative terms

%Calculate u_xx*u_yy, and leave only the nonnegative elements.
u_xy = max(1./h.^2.*(u(1:end-2,2:end-1)+u(3:end,2:end-1)-2.*u(2:end-1,2:end-1)),0)...
       .*max(1./h.^2.*(u(2:end-1,1:end-2)+u(2:end-1,3:end)-2.*u(2:end-1,2:end-1)),0);

if basesNum == 2
    
    %Calculate u_vv*u_v(perp)v(perp), and leave only the nonnegative elements.
    u_vw = max(1./(2.*h.^2).*(u(3:end,1:end-2)+u(1:end-2,3:end)-2.*u(2:end-1,2:end-1)),0)...
           .*max(1./(2.*h.^2).*(u(3:end,3:end)+u(1:end-2,1:end-2)-2.*u(2:end-1,2:end-1)),0);

    %Take the minimum of the two matrices.
    U = min(u_xy,u_vw);
    
else
    
    U = u_xy;
    
end
       
end