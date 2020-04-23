%% 1 load data
clear
clc
close all
%%
numimgs = size(imfinfo('nerves_part.tiff'),1);

for i = 1:numimgs
      Dat(:,:,i) = imread('nerves_part.tiff',i);
      %figure
      %imshow(Dat(:,:,i))
end

%% 2
Nerves = im2double(Dat);
%% 3
start = Nerves(:,:,1);
x0 = 215; y0 = 170;
r = 12+12;
n = 100;
alpha = linspace(0,2*pi,n);

snakex = x0 + r*cos(alpha);
snakey = y0 + r*sin(alpha);
%figure
%imshow(start);
%hold on
%plot(snakex,snakey,'r','linewidth',2)
%plot(snakex,snakey,'ro')
%axis on
%hold off
% 9
C = [snakex' snakey'];

alpha = 0.05*50;
beta = 0.5*10;
Image = start;
N = 30;
tau = 20;
%
slange1 = DerformSegmentationNervNoImg(Image,C,alpha,beta,10,N);
slange1_3d(:,:,1) = slange1;
figure
imshow(Image);
hold on
plot(slange1_3d(:,2,1),slange1_3d(:,1,1),'r','linewidth',2)
hold off
%%
for i = 2:200
    if (i < 6)
        slange1_3d(:,:,i) = DerformSegmentationNervNoImg(Nerves(:,:,i),slange1_3d(:,:,i-1),alpha*1.5,beta,10,10);
    elseif (i < 1000)
        slange1_3d(:,:,i) = DerformSegmentationNervNoImg(Nerves(:,:,i),slange1_3d(:,:,i-1),alpha*1.5,beta,10,10);
    end
end
figure
subplot(1,2,1)
imshow(Image)
hold on
plot(slange1(:,2),slange1(:,1),'b','linewidth',2)
hold off
subplot(1,2,2)
imshow(Nerves(:,:,i));
hold on
plot(slange1_3d(:,2,i),slange1_3d(:,1,i),'r','linewidth',2)
hold off
%%
slange1_3dx = reshape(slange1_3d(:,1,:),100*i,1);
slange1_3dy = reshape(slange1_3d(:,2,:),100*i,1);
for i=1:200
    slange_z((100*(i-1)+1):(100*i)) = ones(1,100)*i;
end
plot3(slange1_3dx,slange1_3dy,slange_z,'b');
%%
x0 = 225; y0 = 95;
r = 12+16;
n = 100;
alpha = linspace(0,2*pi,n);

snakex = x0 + r*cos(alpha);
snakey = y0 + r*sin(alpha);
imshow(start);
hold on
plot(snakex,snakey,'r','linewidth',2)
%plot(snakex,snakey,'ro')
axis on
% 9
C = [snakex' snakey'];

alpha = 1.9;
beta = 1.5;
Image = start;
N = 20;
tau = 20;
%
figure(2)
slange2 = DerformSegmentationNervNoImg(Image,C,alpha,beta,10,N);
slange2_3d(:,:,1) = slange2;
for i = 2:50
    slange2_3d(:,:,i) = DerformSegmentationNervNoImg(Nerves(:,:,i),slange2_3d(:,:,i-1),alpha,beta,10,10);
end
%%
x0 = 115; y0 = 175;
r = 12+14;
n = 100;
alpha = linspace(0,2*pi,n);

snakex = x0 + r*cos(alpha);
snakey = y0 + r*sin(alpha);
imshow(start);
hold on
plot(snakex,snakey,'r','linewidth',2)
%plot(snakex,snakey,'ro')
axis on
%% 9
C = [snakex' snakey'];

alpha = 4.9;
beta = 5.5;
Image = start;
N = 20;
tau = 20;
%
figure(3)
slange3 = DerformSegmentationNervNoImg(Image,C,alpha,beta,10,N);
slange3_3d(:,:,1) = slange3;
for i = 2:20
    slange3_3d(:,:,i) = DerformSegmentationNervNoImg(Nerves(:,:,i),slange3_3d(:,:,i-1),alpha,beta,10,10);
end
%% 2d plot
figure(4)
imshow(start);
hold on
plot(slange1(:,2),slange1(:,1),'r','linewidth',2)
plot(slange2(:,2),slange2(:,1),'g','linewidth',2)
plot(slange3(:,2),slange3(:,1),'b','linewidth',2)
hold off
%%
figure(5)
imshow(Nerves(:,:,20));
hold on
plot(slange1(:,2),slange1(:,1),'r','linewidth',2)
plot(slange1_3d(:,2,20),slange1_3d(:,1,10),'b','linewidth',2);
plot(slange2(:,2),slange2(:,1),'r','linewidth',2)
plot(slange2_3d(:,2,20),slange2_3d(:,1,10),'b','linewidth',2);
plot(slange3(:,2),slange3(:,1),'r','linewidth',2)
plot(slange3_3d(:,2,20),slange3_3d(:,1,10),'b','linewidth',2);
%% 3d plot

slange1_3dx = reshape(slange1_3d(:,1,:),100*20,1);
slange1_3dy = reshape(slange1_3d(:,2,:),100*20,1);
slange2_3dx = reshape(slange2_3d(:,1,:),100*20,1);
slange2_3dy = reshape(slange2_3d(:,2,:),100*20,1);
slange3_3dx = reshape(slange3_3d(:,1,:),100*20,1);
slange3_3dy = reshape(slange3_3d(:,2,:),100*20,1);
for i=1:100
    slange_z((100*(i-1)+1):(100*i)) = ones(1,100)*i;
end
%%
figure(6)
plot3(slange1_3dx,slange1_3dy,slange_z,'b');
hold on
plot3(slange2_3dx,slange2_3dy,slange_z,'r');
plot3(slange3_3dx,slange3_3dy,slange_z,'g');

