%% 1 load data
clear
clc
close all

Dat = VideoReader('data/crawling_amoeba.mov');

i = 1;
while hasFrame(Dat)
    frame = readFrame(Dat);
    fram(:,:,i) = rgb2gray(frame);
    %imshow(fram(:,:,i));
    i = i+1;
end

%% 2
Amoeba = im2double(fram);
%% 3
start = Amoeba(:,:,1);
x0 = 200; y0 = 150;
r = 80;
n = 100;
alpha = linspace(0,2*pi,n);

snakex = x0 + r*cos(alpha);
snakey = y0 + r*sin(alpha);
imshow(start);
hold on
plot(snakex,snakey,'r','linewidth',2)
%plot(snakex,snakey,'ro')
axis on
%% 4
bwIn = poly2mask(snakex,snakey,300,400);
imshow(bwIn);
hold on
plot(snakex,snakey,'r','linewidth',2)

bwOut = ~bwIn;
MeanIn = sum(sum(start.*bwIn))/sum(sum(bwIn));
MeanOut = sum(sum(start.*bwOut))/sum(sum(bwOut));

%% 5
fext = zeros(length(snakex),1);
Inten = zeros(length(snakex),1);
for i = 1:length(snakex)
    fext(i) = (MeanIn - MeanOut)*(2*start(round(snakex(i)),round(snakey(i))) - MeanIn - MeanOut);
    Inten(i) = start(round(snakex(i)),round(snakey(i)));
end
x1 = 1:100;
y1 = ones(100,1)*MeanIn;
y2 = ones(100,1)*MeanOut;
y3 = ones(100,1)*((MeanIn+MeanOut)/2);
plot(x1,Inten,'r','linewidth',2)
hold on;
plot(x1,y1,'b--')
plot(x1,y2,'b--')
plot(x1,y3,'b-')
%% 6
C = [snakex;snakey]';
N = SnakeNormal(C);
imshow(start);
set(gca,'YDir','normal');
hold on
plot(snakex,snakey,'r','linewidth',2)
tau = 100;
for i = 1:100
    plot([snakex(i),snakex(i)+N(i,1)*fext(i)*tau],[snakey(i),snakey(i)+N(i,2)*fext(i)*tau],'b-','linewidth',2) 
end
axis on

%% 7
n = 100;
alpha = 0.5;
beta = 0.5;
Bint = ImplicitSmoothMat(alpha,beta,n);

tau = 20;
con = tau*diag(fext)*N;
snakenew = [snakex' + con(:,1),snakey' + con(:,2)];

snak = Bint*snakenew;
imshow(start);
hold on
%plot(snakenew(:,1),snakenew(:,2),'r','linewidth',2)
%plot(snakex,snakey,'b','linewidth',2)
plot(snak(:,1),snak(:,2),'r','linewidth',2)

%% 8 

%% 9
x0 = 200; y0 = 150;
r = 80;
n = 100;
alpha = linspace(0,2*pi,n);

snakex = x0 + r*cos(alpha);
snakey = y0 + r*sin(alpha);

C = [snakex' snakey'];

alpha = 0.05;
beta = 0.5;
Image = start;
N = 100;
tau = 20;
slange = DerformSegmentation(Image,C,alpha,beta,tau,N);





