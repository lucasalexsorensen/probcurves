function SmotMat = ImplicitSmoothMat(alpha,beta,n)

I = eye(n);
A = zeros(n,n);
B = zeros(n,n);
A(1,1) = -2;
A(1,2) = 1;
A(2,2) = -2;
A(2,2+1) = 1;
A(2,2-1) = 1;
B(1,1) = -6;
B(1,2) = 4;
B(1,3) = -1;
B(2,1) = 4;
B(2,2) = -6;
B(2,3) = 4;
B(2,4) = -1;

for i = 3:(n-2)
    A(i,i) = -2;
    A(i,i+1) = 1;
    A(i,i-1) = 1;
    
    B(i,i) = -6;
    B(i,i+1) = 4;
    B(i,i+2) = -1;
    B(i,i-1) = 4;
    B(i,i-2) = -1;
end
A(n-1,n-1) = -2;
A(n-1,n) = 1;
A(n,n-1) = 1;
A(n,n) = -2;
A(n-1,n-2) = 1;

B(n,n) = -6;
B(n,n-1) = 4;
B(n,n-2) = -1;
B(n-1,n) = 4;
B(n-1,n-1) = -6;
B(n-1,n-2) = 4;
B(n-1,n-3) = -1;
%Corner
A(1,n) = 1;
A(n,1) = 1;
B(1,n) = 4;
B(n,1) = 4;
B(1,n-1) = -1;
B(n,2) = -1;
B(2,n) = -1;
B(n-1,1) = -1;

SmotMat = inv(I-(alpha*A+beta*B));