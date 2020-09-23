function [estimated_variance, residuals] = var_est_DCT2D(input_image,blk_size,degre,mode)
% function [ estimated_variance , estimated_expectation , residuals ] = var_est_DCT2D(input_image,blk_size,degre,mode)
% function [ estimated_variance , estimated_expectation ] = var_est_DCT2D(input_image,blk_size,degre,mode)
% Estimation of expectation and variance of pixels based on a 2D-DCT
% (trigonometric polynomial) model
% ----------------------------------------------------------------------------------------------- %
%		Input Parameters :
%	input_image			Matrix : image to be inspected
%	blk_size			one or two values : specify the size of blocks
%	degree				int pos: number of coefficients used
%	mode				1/2    : specify whether variance is estimated per
%                                block or for each pixel
% -----------------%
%		Output parameters :
%	estimated_variance		matrix : estimation of variance
%   estimated_expectation	matrix : estimation of expectation
%   residuals				matrix : residuals (input_img-est_expectation)
% ----------------------------------------------------------------------------------------------- %
% Note that the code is provided 'as is' without any warranty,
% Implementation in design for grayscale image, a loop on the thrid dimension should do the extension for color images
%
% Programmer : Rémi Cogranne
% Written for Matlab 7.14.0.739 (R2012a)
%

% Verification and integrity of input arguments 
if nargin < 1;
	error('Input argument "input_image", that is the input image, is undefined.');
end,

if nargin < 2;
	blk_size = 4;
end,

if nargin < 3;
	degre = 6;
end,

if nargin < 4,
    mode = 1;
end,

if numel(blk_size) == 1,
    blk_size_x = blk_size;
    blk_size_y = blk_size;
elseif numel(blk_size) == 2
    blk_size_x = blk_size(1);
    blk_size_y = blk_size(2);
else    
   	error('Input argument "block_size", must contains one or two value only !');
end

input_image = double(input_image);
% Resize input_image to insure it contains Kn blocks

q = (degre+1)*(degre+2)/2;                % q = number of parameters per block = sum of i for i=0:degre+1
[Mx My] = size(input_image);

if degre>blk_size_y || degre>blk_size_x, error('Number of basis vectors excesses block dimension !!'); end

tmp = zeros(blk_size_x , blk_size_y);
k=1;
for kx = 0:degre;
    for ky = 0:degre-kx;
        if kx < blk_size_x && ky < blk_size_y
            tmp2 = tmp; tmp2(kx+1,ky+1) = 1;
            vec_tmp = idct2(tmp2);
            H(:,k) = vec_tmp(:);
            k=k+1;
        end,
    end,
end

% WF = fspecial('gaussian',blk_size,Sigma);
% W = diag(WF(:));
% PHorth = eye(blk_size_x*blk_size_y) - H*inv(H'*W*H)*H'*W;
% PH = H*inv(H'*W*H)*H'*W;

PHorth = eye(blk_size_x*blk_size_y) - H*inv(H'*H)*H';
PH = H*inv(H'*H)*H';


% if mode==1 then variance is estimated block wise
if mode==1,

    % Resize input_image to insure it contains Kn blocks
    input_image = input_image(1:end-mod(size(input_image,1),blk_size_x),1:end-mod(size(input_image,1),blk_size_y) );
    estimated_variance = zeros(size(input_image));       % Pixel variance
    estimated_expectation = zeros(size(input_image));       % Pixel variance
    Kn = floor(Mx/blk_size_x) * floor(My/blk_size_y);          % Number of blocks of m pixels
    Y = zeros(blk_size_x * blk_size_y,Kn);          % Y will hold block pixel values as columns

    %Here, the input image is first reshaped
    k=1;
    for i = 1 : blk_size_x, 
        for j = 1 : blk_size_y, % Form Kn blocks of blk_size_x*blk_size_y pixels as
            aux = input_image(i:blk_size_x:end,j:blk_size_y:end);   % columns of Y
            Y(k,:) = aux(:)';
            k=k+1;
        end
    end

    % then the variance and expectation are calculated for each subblock
    sig2 = PHorth * Y;
    sig2 = sum(sig2.^2) / (blk_size_x*blk_size_y-q);% variance in kth block

    tmp_expectation = PH * Y;

    Sy = ones(blk_size_x*blk_size_y,1) * sig2;                        % Variance of all pixels (order as in Y)

    k=1;
    for i = 1 : blk_size_x             % Reshaping the variance Sy to size of input_image
        for j = 1 : blk_size_y,
            estimated_variance(i:blk_size_x:end,j:blk_size_y:end) = reshape(Sy(k,:),size(input_image(i:blk_size_x:end,j:blk_size_y:end)));
            estimated_expectation(i:blk_size_x:end,j:blk_size_y:end) = reshape(tmp_expectation(k,:),size(input_image(i:blk_size_x:end,j:blk_size_y:end)));
        end,
    end

end,

% if mode==2 then variance is estimated for each pixel
if mode==2,

    estimated_variance = zeros(size(input_image));       % Pixel variance
    estimated_expectation = zeros(size(input_image));       % Pixel variance

    if ( mod(blk_size_x,2) ~= 1 || mod(blk_size_y,2)~= 1 ) ; error('The block dimension should be odd numbers (or you should modify the code)'); end
    hx = ceil( (blk_size_x-1)/2);
    hy = ceil( (blk_size_y-1)/2);
    c = ceil( (blk_size_x*blk_size_y+1)/2);


    for i=hx+1:size(input_image,1)-hx-1;
        for j=hy+1:size(input_image,1)-hy-1;
            tmp = input_image(i-hx:i+hx , j-hy:j+hy);
            
            estimated_variance(i,j) = norm( PHorth * tmp(:) )^2/(blk_size_x * blk_size_y - q);
            %The two following lines can be removed for speed (hence, the
            %estimated expectation is not given)
            tmp = PH * tmp(:);
            estimated_expectation(i,j) = tmp(c);
        end,
    end

end,

residuals = estimated_expectation - input_image;

return,
