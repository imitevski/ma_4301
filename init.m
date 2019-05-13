function [u] = init(F,g,n,h,X,Y)
%init.m calculates a reasonable initial guess to plug into the FAS
%function.

m = n-2;

F = sqrt(2*F(2:m+1,2:m+1));

G = g(X,Y);

F(:,1) = F(:,1) - G(2:m+1,1)/h^2;
F(:,m) = F(:,m) - G(2:m+1,m+2)/h^2;
F(1,:) = F(1,:) - G(1,2:m+1)/h^2;
F(m,:) = F(m,:) - G(m+2,2:m+1)/h^2;

F = reshape(F,m^2,1);

I = speye(m);
e = ones(m,1);
T = spdiags([e -4*e e],[-1 0 1],m,m);
S = spdiags([e e],[-1 1],m,m);
A = (kron(I,T) + kron(S,I))/h^2;

uVec = A\F;
G(2:m+1,2:m+1) = reshape(uVec,[m,m]);

u = G;

end