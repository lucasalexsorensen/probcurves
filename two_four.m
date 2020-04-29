clear all;
addpath('probabilistic_data');

% load image
im = imread('simple_test.png');
im_cw = reshape(im, [size(im,1)*size(im,2), 1]); % columnwise

% init snake
[xs, ys] = build_snake(250, 250, 175);
mask = poly2mask(xs, ys, size(im, 1), size(im, 2));
c_in = reshape(mask, [size(im,1)*size(im,2), 1]);
A_in = length(find(mask > 0));
A_out = size(im,1)*size(im,2);

% build biadjacency matrix
B = build_B(im_cw);

% compute freqs
f_in = (B'*c_in)/A_in;
f_out = (B'*~c_in)/A_out;

function B = build_B (im_cw)
  B = zeros(length(im_cw), 256);
  for i=1:length(im_cw)
    B(i, im_cw(i) + 1) = 1;
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