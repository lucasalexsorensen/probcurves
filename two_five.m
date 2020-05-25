clear all;
addpath('probabilistic_data');

% load image
im = imread('textured_test.png');
im = imgaussfilt(im, 1);
im_cw = reshape(im, [size(im,1)*size(im,2), 1]); % columnwise unwrap

% make patches
M = 11; % patch size
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
bins = 30; % amount of clusters
[idxs, centroids] = kmeans(vectors, bins);
% centroids = round(centroids); % round here? since pixel values are not doubles

% gather centroids as a long column vector of patches
dict_patches = reshape(centroids, [M*M*bins, 1]); % 

%%
% build assignment image S

im = 255.*im2double(im);
pad = ceil(M/2);
S = zeros(size(im));
tot = size(S,1)*size(S,2);
for i=1:size(S,1)
  for j=1:size(S,2)
    if (mod(j, 500) == 0)
      fprintf('%f percent done\n', 100*((i-1)*size(S,1) + j) / tot)
    end
    if i < pad || j < pad || i > size(S,1) - pad || j > size(S,2) - pad
      S(i,j) = 0;
      continue
    end
    
    patch = reshape(im(i-(pad-1):i+(pad-1), j-(pad-1):j+(pad-1)), [1, M*M]);
    S(i,j) = knnsearch(centroids, patch);
  end
end

imagesc(S);


%%
% build biadjacency matrix B
I = [];
J = [];
V = [];
S_cw = reshape(S, [size(S,1)*size(S,2), 1]); % columnwise unwrap
for i=1:size(S,1) * size(S,2)
  for j=1:length(dict_patches)
    if S_cw(i) == j
      I(end+1) = i; % row indices
      J(end+1) = j; % col indices
      V(end+1) = 1; % values
    end
  end
end

B = sparse(I, J, V, size(im,1) * size(im,2), length(dict_patches));
global rowSums;
rowSums = sum(B, 2);

%%
% build snake
[xs, ys] = build_snake(330, 330, 120);
C = [xs' ys'];

mask = poly2mask(xs, ys, size(im, 1), size(im, 2)); % discretizes snake to mask im
c_in = reshape(mask, [size(im,1)*size(im,2), 1]); % reshapes mask as a column
A_in = length(find(mask > 0)); % area of mask
A_out = size(im,1)*size(im,2) - A_in; % area of non-mask (total area - mask area)
f_in = (B'*c_in); % frequency of pixels in mask
f_out = (B'*~c_in); % frequency of pixels outside mask
Z = f_in + f_out; % normalization const
p_in = f_in./Z; % normalize such that probs are between 0-1
p_in(isnan(p_in)) = 0.5; % if nan (e.g. both freqs are 0), set to 0.5

%% visualize patches + patch probs
subplot(1,2,1)
patchvis = reshape(dict_patches, [bins, M*M]);
patchim = zeros(M*sqrt(bins), M*sqrt(bins));
for i=1:sqrt(bins)
  for j=1:sqrt(bins)
    row = 1 + (i-1) * M;
    col = 1 + (j-1) * M;
    p = reshape(patchvis((i-1) * sqrt(bins) + j, :), [M M]);
    patchim(row:(row+M-1), col:(col+M-1)) = p;
  end
end
imagesc(patchim);
colormap gray;

subplot(1,2,2)
patchprob = reshape(p_in, [bins, M*M]);
colors = mean(patchprob, 2);
colorim = zeros(M*sqrt(bins), M*sqrt(bins));
colors_ij = reshape(colors, [sqrt(bins) sqrt(bins)]);
for i=1:sqrt(bins)
  for j=1:sqrt(bins)
    row = 1 + (i-1) * M;
    col = 1 + (j-1) * M;
    colorim(row:(row+M-1), col:(col+M-1)) = ones(M, M) * colors_ij(i, j);
  end
end
imagesc(colorim);

%%
P_in = reshape((B*p_in) ./ rowSums, [size(im,1), size(im,2)]);
imagesc(P_in);
colormap redblue;
hold on;
plot([C(:,2); C(1,2)],[C(:,1); C(1,1)],'r','linewidth',2)


%%
[xs, ys] = build_snake(140, 330, 120);
C = [xs' ys'];


filename = 'testPatchProbs.gif';
smoothMat = ImplicitSmoothMat(0.0015, 0.0001, length(xs));
for i = 1:200
  tau = 10;
  
  if (i == 40)
    tau = 2;
  end
  
  if (i == 75)
    tau = tau / 8;
  end
  [C, P_in_iteration] = iterate(im, smoothMat, B, [xs' ys'], tau);
  xs = C(:,1)';
  ys = C(:,2)';
  if mod(i, 5) == 0
    h = figure();
    subplot(1,2,1);
    imagesc(P_in_iteration);
    axis image;
    
    subplot(1,2,2);
    imshow(im./255);
    hold on;
    plot([C(:,2); C(1,2)],[C(:,1); C(1,1)],'r','linewidth',2)
    title(sprintf('iteration %d', i));

    frame = getframe(h);
    im2 = frame2im(frame);
    [imind,cm] = rgb2ind(im2,256); 
    if i == 5 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end
  end
end


function [C, P_in] = iterate (im, smoothMat, B, C, tau)
  global rowSums
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
  P_in = reshape((B*p_in) ./ rowSums, [size(im, 1), size(im, 2)]);
  P_in(isnan(P_in)) = 0.5;
  P_out = 1-P_in;

  fext = zeros(length(xs),1);
  for j = 1:length(xs)
    fext(j) = P_in(round(xs(j)),round(ys(j))) - P_out(round(xs(j)),round(ys(j)));
  end
  normals = SnakeNormal([xs' ys']);

  con = tau*diag(fext)*normals;
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
