clear all;
addpath('probabilistic_data');

% load image
global imsize masksize;
im = imread('simple_test.png');
imsize = (size(im, 1) * size(im, 2));

% init snake
[xs, ys] = build_snake(250, 250, 175);
mask = poly2mask(xs, ys, size(im, 1), size(im, 2));
masksize = length(find(im(mask)));
x
global probcache;
probcache = -1 * ones(1, 256);

for i = 1:size(im, 1)
  for j = 1:size(im, 2)
    p_in = probs(im, mask, im(i,j));
    im(i,j) = p_in * 255;
  end
end

imshow(im);

function [p_in, p_out] = probs (im, mask, intensity)
  global imsize masksize probcache;
  if probcache(intensity+1) ~= -1
    p_in = probcache(intensity+1);
    p_out = 1-p_in;
  else
    f_in = length(find(im(mask) == intensity)) / masksize;
    f_out = length(find(im(~mask) == intensity)) / imsize;
    norm_const = 1 / (f_in + f_out);
    p_in = norm_const * f_in;
    p_out = norm_const * f_out;
    probcache(intensity+1) = p_in;
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