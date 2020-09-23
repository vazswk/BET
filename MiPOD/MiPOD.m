function cost = MiPOD( Cover,alpha )

% Reading cover image
% Cover = double(imread('1.pgm'));

% Set MVG parametrs:
% -  Payload 0.2 bpp
% -  Variance estimator: DCT2D-9-8
% -  Pre-filter cover with wiener filter on 2*2 support
% -  Post-filter Fisher info with average filtering on 7*7 support
EmbeddingParam.Payload = alpha;
EmbeddingParam.Filter.ImagePreFilter = true;
EmbeddingParam.Filter.ImagePreFilterType = 'Wiener';
EmbeddingParam.Filter.ImagePreFilterC1 = 2;
EmbeddingParam.Filter.ImagePreFilterC2 = 2;
EmbeddingParam.Filter.VarPostFilter = false;
EmbeddingParam.Filter.CostsPostFilter = false;
EmbeddingParam.Filter.FIPostFilter = true;
EmbeddingParam.Filter.FIPostFilterType = 'average';
EmbeddingParam.Filter.FIPostFilterSize = 5;
EmbeddingParam.EstimatorSet = {{'DCT2D',7,6}};

% Embed in cover image
% tStart = tic;
% [Stego, ChangeRate, pChange, Cost] = emb_minKL_ternary_multivar_spatial(Cover, EmbeddingParam);
[~, ~, ~, cost] = emb_minKL_ternary_multivar_spatial(Cover, EmbeddingParam);
% tEnd = toc(tStart);

% fprintf('Time to finish embedding: %.4f seconds\n', tEnd);
% 
% % Plot the results
% figure;imshow(Cover-Stego,[]);title('Changed Pixels');
% figure;imshow(pChange,[]);title('Probability of Change for Each Pixel');
% figure;imshow(Cost,[]);title('Embedding Cost at Each Pixel');
end