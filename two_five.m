clear all;
addpath('probabilistic_data');

% load image
im = imread('textured_test.png');
im = imresize(im, 0.4);
im_cw = reshape(im, [size(im,1)*size(im,2), 1]); % columnwise unwrap

% make patches
M = 16; % patch size
patches = extractPatches(im, [M M]);
n_patches = min(floor(size(im,1) / M), floor(size(im,2) / M));

% collect patch vectors
vectors = zeros(n_patches*n_patches, M*M);
for i=1:n_patches
  for j=1:n_patches
    patch = patches{i, j};
    vectors((i-1)*n_patches + j, :) = reshape(patch, [M*M, 1]); %columnwise unwrap
  end
end

% k-means clustering on patch vectors
bins = 50; % amount of clusters
[idxs, centroids] = kmeans(vectors, bins);
centroids = round(centroids); % round here? since pixel values are not doubles

% gather centroids as a long column vector of patches
dict_patches = reshape(centroids, [M*M*bins, 1]); % 

% build biadjacency matrix B
B = zeros(size(im,1) * size(im,2), length(dict_patches), 'single');
for i=1:size(B,1)
  for j=1:size(B,2)
    if im_cw(i) == dict_patches(j)
      B(i,j) = 1;
    end
  end
end
rowSums = sum(B, 2);
%%
% build snake
[xs, ys] = build_snake(130, 110, 50);

%P_in = (B*p_in)./rowSums;
%imagesc(reshape(P_in, [size(im,1), size(im,2)]));

imshow(im);
hold on;
smoothMat = ImplicitSmoothMat(0.01, 0.3, length(xs));
for i = 1:10
  C = iterate(im, smoothMat, B, [xs' ys']);
  xs = C(:,1)';
  ys = C(:,2)';
  plot([C(:,2); C(1,2)],[C(:,1); C(1,1)],'r','linewidth',2)
  drawnow;
  disp(i);
end


%%
function C = iterate (im, smoothMat, B, C)
  xs = C(:,1)';
  ys = C(:,2)';

  mask = poly2mask(xs, ys, size(im, 1), size(im, 2));
  c_in = reshape(mask, [size(im,1)*size(im,2), 1]);
  A_in = length(find(mask > 0));
  A_out = size(im,1)*size(im,2) - A_in;
  f_in = (B'*c_in)/A_in;
  f_out = (B'*~c_in)/A_out;
  Z = f_in + f_out;
  p_in = f_in./Z;
  p_in(isnan(p_in)) = 0.5;
  p_out = 1 - p_in;

  fext = zeros(length(xs),1);
  for j = 1:length(xs)
      snakeim = im(round(xs(j)),round(ys(j)));
      fext(j) = p_in(snakeim+1) - p_out(snakeim+1);
  end
  normals = SnakeNormal([xs' ys']);

  con = 10*diag(fext)*normals;
  snakenew = [xs' + con(:,1), ys' + con(:,2)];

  C = smoothMat*snakenew;
  S = distribute_points(C);
  C = remove_intersections(S);
end


function patches = extractPatches (im, patchSz)
  imSz = size(im);
  xIdxs = [1:patchSz(2):imSz(2) imSz(2)+1];
  yIdxs = [1:patchSz(1):imSz(1) imSz(1)+1];
  patches = cell(length(yIdxs)-1,length(xIdxs)-1);
  for i = 1:length(yIdxs)-1
      Isub = im(yIdxs(i):yIdxs(i+1)-1,:);
      for j = 1:length(xIdxs)-1
          patches{i,j} = Isub(:,xIdxs(j):xIdxs(j+1)-1);
      end
  end
end

function [xs, ys] = build_snake (x0, y0, r, n)
  if nargin < 4
    n = 50;
  end
  alpha = linspace(0, 2*pi, n);
  xs = x0 + r*cos(alpha);
  ys = y0 + r*sin(alpha);
end

function reconst = reconstFromVectors (vectors, imSize, n_patches, M)
  reconst = zeros(imSize, imSize);
  for i=1:n_patches
    for j=1:n_patches
      row = 1 + (i-1) * M;
      col = 1 + (j-1) * M;
      reconst(row:(row+M-1), col:(col+M-1)) = reshape(vectors((i-1)*n_patches + j, :), [M M]);
    end
  end
end
