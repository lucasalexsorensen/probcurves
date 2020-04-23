function [slange,PimStart,PimEnd] = DerformSegmentation(Image,C,alpha,beta,tau,N1,lim)
if(lim>N1)
    lim = N1;
end
im = Image;
im = round(im);
n = size(C,1);
[Nim,Mim] = size(im);
Bint = ImplicitSmoothMat(alpha,beta,n);
%figure(1)
imshow(im./255);
hold on
plot(C(:,2),C(:,1),'r','linewidth',1)
hold off
imVec = reshape(im,Nim*Mim,1);
B = zeros(length(imVec),256);
for k = 1:256
    B(imVec==(k-1),k) = 1;
end


for i = 1:N1
    if(i == lim)
        tau = 5;
    end
    snakex = C(:,1);
    snakey = C(:,2);
    
    % Step 4    
    bwIn = poly2mask(snakex,snakey,Nim,Mim);
    bwOut = ~bwIn;
    bwIn = ~bwOut;
    cIn = reshape(bwIn,Nim*Mim,1);
    cOut = 1 - cIn;
    AIn = sum(cIn);
    AOut = sum(cOut);
    
    fin = B'*cIn./AIn;
    fout = B'*cOut./AOut;
    
    
    Z = fin + fout;
    pin = fin./Z;
    pin(isnan(pin)) = 0.5;
    pout = fout./Z;
    pout(isnan(pout)) = 0.5;
    %figure(2)
    %plot(0:255,pin)
    %xlim([0 255])
    
    
    if(i==1)
        PimStart = zeros(size(im));
        PoumStart = zeros(size(im));
        for k=1:256
            idx = im==k;
            PimStart(idx) = pin(k);
            PoumStart(idx) = pout(k);
        end
    elseif(i==N1)
        PimEnd = zeros(size(im));
        PoumEnd = zeros(size(im));
        for k=1:256
            idx = im==k;
            PimEnd(idx) = pin(k);
            PoumEnd(idx) = pout(k);
        end
    end
    % Step 5
    fext = zeros(length(snakex),1);
    for j = 1:length(snakex)
        snakeim = im(round(snakex(j)),round(snakey(j)));
        fext(j) = pin(snakeim+1) - pout(snakeim+1);
    end

    % Step 6
    N = SnakeNormal(C);
    
    % Step 7
    con = tau*diag(fext)*N;
    snakenew = [snakex + con(:,1),snakey + con(:,2)];

    C = Bint*snakenew;
    % Step 8
    S = distribute_points(C);
    C = remove_intersections(S);
    % Plotting
    figure(1)
    imshow(im./255);
    title("n = " + i)
    hold on
    plot([C(:,2); C(1,2)],[C(:,1); C(1,1)],'r','linewidth',2)
    %for k = 1:length(snakex)
    %    plot([C(k,2),C(k,2)+N(k,2)*fext(k)*tau],[C(k,1),C(k,1)+N(k,1)*fext(k)*tau],'b-','linewidth',2) 
    %end
    %disp(i);
end
% Output
slange = C;