function N = SnakeNormal(C)

n = size(C,1);
N = zeros(size(C));

c1 = C(2,:)-C(n,:);
N(1,:) = -[-c1(2) c1(1)]./norm(c1,2);
for i = 2:(n-1)
    c1 = C(i+1,:)-C(i-1,:);
    N(i,:) = -[-c1(2) c1(1)]./norm(c1,2);
end
c1 = C(n-1,:)-C(1,:);
N(n,:) = -[-c1(2) c1(1)]./norm(c1,2);