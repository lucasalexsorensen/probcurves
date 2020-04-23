%% 1 load data
clear
clc
close all

Dat = imread('Data/simple_test.png');



%% 2
start = double(Dat);
%% 3
close all
x0 = 210; y0 = 210;
r = 170;
n = 200;
alpha = linspace(0,2*pi,n);

snakex = x0 + r*cos(alpha);
snakey = y0 + r*sin(alpha);
imshow(start./256);
hold on
plot(snakex,snakey,'r','linewidth',2)
%plot(snakex,snakey,'ro')
axis on

C = [snakex' snakey'];

alpha = 0.01;
beta = 0.3;
Image = start;
N = 40;
tau = 20;
[slange,PinStart,PinEnd] = DerformSegmentation(Image,C,alpha,beta,tau,N);
%%
figure
subplot(1,2,1)
imagesc(PinStart)
subplot(1,2,2)
imagesc(PinEnd)
%% New image
Dat = imread('Data/overlap_test.png');
overlap = double(Dat);
close all
x0 = 210; y0 = 210;
r = 170;
n = 200;
alpha = linspace(0,2*pi,n);

snakex = x0 + r*cos(alpha);
snakey = y0 + r*sin(alpha);
imshow(overlap./256);
hold on
plot(snakex,snakey,'r','linewidth',2)
%plot(snakex,snakey,'ro')
axis on

%%
C = [snakex' snakey'];

alpha = 0.01;
beta = 0.3;
N = 90;
tau = 20;
lim = 70;
[slange,PinStart,PinEnd] = DerformSegmentation(overlap,C,alpha,beta,tau,N,lim);

%%
figure
subplot(1,2,1)
imagesc(PinStart)
subplot(1,2,2)
imagesc(PinEnd)
%% Image 3
Dat = imread('Data/textured_test.png');
textured = double(Dat);
close all
x0 = 210; y0 = 210;
r = 170;
n = 200;
alpha = linspace(0,2*pi,n);

snakex = x0 + r*cos(alpha);
snakey = y0 + r*sin(alpha);
imshow(textured./256);
hold on
plot(snakex,snakey,'r','linewidth',2)
%plot(snakex,snakey,'ro')
axis on

%%
C = [snakex' snakey'];

alpha = 0.01;
beta = 0.3;
N = 90;
tau = 20;
lim = 70;
[slange,PinStart,PinEnd] = DerformSegmentation(textured,C,alpha,beta,tau,N,lim);

%%
figure
subplot(1,2,1)
imagesc(PinStart)
subplot(1,2,2)
imagesc(PinEnd)


