function slange = DerformSegmentation(Image,C,alpha,beta,tau,N1)
im = Image;
n = size(C,1);

Bint = ImplicitSmoothMat(alpha,beta,n);
imshow(im);
hold on
plot(C(:,2),C(:,1),'r','linewidth',1)
for i = 1:N1
    snakex = C(:,1);
    snakey = C(:,2);
    
    % Step 4
    bwIn = poly2mask(snakex,snakey,300,400);
    bwOut = ~bwIn;
    MeanIn = sum(sum(im.*bwIn))/sum(sum(bwIn));
    MeanOut = sum(sum(im.*bwOut))/sum(sum(bwOut));
    
    % Step 5
    fext = zeros(length(snakex),1);
    for j = 1:length(snakex)
        fext(j) = (MeanIn - MeanOut)*(2*im(round(snakex(j)),round(snakey(j))) - MeanIn - MeanOut);
    end

    % Step 6
    N = SnakeNormal(C);
    
    % Step 7
    con = tau*diag(fext)*N;
    snakenew = [snakex + con(:,1),snakey + con(:,2)];

    C = Bint*snakenew;
    % Plotting
    figure(1)
    imshow(im);
    plot(C(:,2),C(:,1),'r','linewidth',2)
    for k = 1:100
        plot([snakey(k),snakey(k)+N(k,2)*fext(k)*tau],[snakex(k),snakex(k)+N(k,1)*fext(k)*tau],'b-','linewidth',2) 
    end
    
    % Step 8
    S = distribute_points(C);
    C = remove_intersections(S);
end
% Output
slange = C;