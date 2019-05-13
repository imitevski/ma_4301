function [u] = notJacobi(u,F,h,basesNum)
% notJacobi performs Gauss-Seidel iterations

indices = zeros(size(u));
indices(2:end-1,2:end-1) = 1;

indices = find(indices == 1);
    
for k = 1:length(indices)

    [i,j] = ind2sub(size(u),indices(k));
          
    if basesNum == 2
        
        A_xy = (1/h)^4*(u(i-1,j)+u(i+1,j)-2*u(i,j))...
                      *(u(i,j-1)+u(i,j+1)-2*u(i,j));
          
        A_vw = 1/(4*h^4)*(u(i-1,j-1)+u(i+1,j+1)-2*u(i,j))...
                        *(u(i+1,j-1)+u(i-1,j+1)-2*u(i,j));

        if A_xy <= A_vw

            u(i,j) = 0.25*(u(i+1,j)+u(i-1,j)+u(i,j-1)+u(i,j+1))...
                -0.5*sqrt(0.25*((u(i+1,j)+u(i-1,j)-u(i,j-1)-u(i,j+1))^2)...
                +h^4*F(i,j));

        else

            u(i,j) = 0.25*(u(i-1,j+1)+u(i+1,j-1)+u(i+1,j+1)+u(i-1,j-1))...
                -0.5*sqrt(0.25*((u(i-1,j+1)+u(i+1,j-1)-u(i+1,j+1)-u(i-1,j-1))^2)...
                +4*h^4*F(i,j));

        end
        
    else
        
            u(i,j) = 0.25*(u(i+1,j)+u(i-1,j)+u(i,j-1)+u(i,j+1))...
                -0.5*sqrt(0.25*((u(i+1,j)+u(i-1,j)-u(i,j-1)-u(i,j+1))^2)...
                +h^4*F(i,j));
        
    end

end

u = u(2:end-1,2:end-1);

end