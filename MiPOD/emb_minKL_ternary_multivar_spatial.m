function [Y, change_rate, pChange, Cost] = emb_minKL_ternary_multivar_spatial(X,EmbeddingParam)
%
% Inputs :
% X - is the input cover image.
% EmbeddingParam -  is a structure which contains all of the parametrs in
% regrads with embedding. The fields are as follows:
% EmbeddingParam
%   .Payload : It is the required payload to be embed in the image in bit 
%              per pixel (bpp)  
%   .EstimatorSet : It contains the set of variance estimators of MVG model
%                   which is a structure that contains a structure with 
%                   three inputs:
%                   - "Estimator Type" : {'POLY2D','DCT2D','LocalVar'}
%                   - "Block Size" : any odd natural number 
%                   - "Degree" or Type" : for any of the first two types 
%                     of estimators it can be any integer from 1 to Block 
%                     Size and for the last estimator it sets the type of
%                     directional filters form {1, 2, ...,20}
%   .Filter  : This field is compounded of several other fileds which can 
%              be used to do the filtering Process during different stages
%              of embedding. It has the following fields:
%      .ImagePreFilter : Used to pre-filter the input cover image before 
%                        estimeting the varicance. Possible inputs are 
%                        "true" or "false" based on whether it is required 
%                        to pre-filter the image
%      .ImagePreFilterType : Possible filters are
%                            {'Wiener','Wavelet','KB','KV'}. For using two
%                            first filters block size dimensions should be
%                            set using two separate following fields:
%                            .ImagePreFilterC1 and .ImagePreFilterC2 which
%                            can be any natural number 
%      .VarPostFilter : Used to post-filter estiated variances. Possible  
%                       inputs are "true" or "false"
%      .VarPostFilterType : Possible filters are 
%                           {'gaussian','average','anisodiff'}. For using
%                           any of the two first filters block size should
%                           be set using the field .VarPostFilterSize which
%                           can be any natural number and for using
%                           Gaussian filter standard diviation of the
%                           filter also needs to be set using the field
%                           .VarPostFilterSigma which can be any real number  
%      .FIPostFilter : Used to post-filter Fisher information field. Possible  
%                       inputs are "true" or "false"
%      .FIPostFilterType : Possible filters are 
%                           {'gaussian','average','anisodiff'}. For using
%                           any of the two first filters block size should
%                           be set using the field .FIPostFilterSize which
%                           can be any natural number and for using
%                           Gaussian filter standard diviation of the
%                           filter also needs to be set using the field
%                           .FIPostFilterSigma which can be any real number  
%      .CostsPostFilter : Used to post-filter computed costs. Possible  
%                       inputs are "true" or "false"
%      .CostsPostFilterType : Possible filters are 
%                           {'gaussian','average','anisodiff'}. For using
%                           any of the two first filters block size should
%                           be set using the field .CostsPostFilterSize which
%                           can be any natural number and for using
%                           Gaussian filter standard diviation of the
%                           filter also needs to be set using the field
%                           .CostsPostFilterSigma which can be any real number
% Outputs:
% Y - Output stego image
% change_rate - Change rate for computing the output stego image
% pChange - Probability of changing each pixel in order to compute the output 
%           stego image which carries certian payload set by EmbeddingParam.Payload
% Cost - Computed costs of embedding at each pixel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ischar( X )
    X = double(imread(X));
else
    X = double( X );
end

alpha = EmbeddingParam.Payload * log(2);
variance_low_limit = 0.01;

if EmbeddingParam.Filter.ImagePreFilter
    switch (EmbeddingParam.Filter.ImagePreFilterType)
        case 'Wiener'
            XRes = {X-wiener2(X,[EmbeddingParam.Filter.ImagePreFilterC1 EmbeddingParam.Filter.ImagePreFilterC2])};
        case 'Wavelet'
            XRes = {X-wavden2(X,EmbeddingParam.Filter.ImagePreFilterC1,EmbeddingParam.Filter.ImagePreFilterC2)};
        case 'KB'
            F_KB = [-0.25, 0.5, -0.25; 0.5, -1, 0.5; -0.25, 0.5, -0.25];
            R = conv2(X, F_KB, 'same');
            XRes = {R};
        case 'KV'
            F_KV = (1/12).*[-1, 2, -2, 2, -1; 2, -6, 8, -6, 2; -2, 8, -12, 8, -2; 2, -6, 8, -6, 2; -1, 2, -2, 2, -1];
            R = conv2(X, F_KV, 'same');
            XRes = {R};
    end            
else
    XRes = {X};
end
%XRes = {NoiseExtractFromImage(X,NoiseSigma)};
FisherInformation = 0;
for R = 1 : numel(XRes)
    for E = 1 : numel(EmbeddingParam.EstimatorSet)
        PadSize = ceil(EmbeddingParam.EstimatorSet{E}{2}/2);
        PaddedXRes = padarray(XRes{R},[PadSize PadSize],'symmetric');
        switch EmbeddingParam.EstimatorSet{E}{1}
            case 'LocalVar'
                multi_var_est = localvar_adaptive(PaddedXRes,EmbeddingParam.EstimatorSet{E}{2},EmbeddingParam.EstimatorSet{E}{3});
            case 'DCT2D'
                multi_var_est = var_est_DCT2D(PaddedXRes,EmbeddingParam.EstimatorSet{E}{2},EmbeddingParam.EstimatorSet{E}{3},2);
            case 'POLY2D'
                multi_var_est = var_est_poly2D(PaddedXRes,EmbeddingParam.EstimatorSet{E}{2},EmbeddingParam.EstimatorSet{E}{3},2);
            otherwise
                multi_var_est = padarray(EmbeddingParam.InputVariance,[PadSize PadSize],'symmetric');
        end
        multi_var_est = multi_var_est(PadSize+1:end-PadSize,PadSize+1:end-PadSize);
        multi_var_est( multi_var_est < variance_low_limit ) = variance_low_limit;
        if EmbeddingParam.Filter.VarPostFilter
            switch (EmbeddingParam.Filter.VarPostFilterType)
                case 'gaussian'
                    h = fspecial(EmbeddingParam.Filter.VarPostFilterType,EmbeddingParam.Filter.VarPostFilterSize,EmbeddingParam.Filter.VarPostFilterSigma);
                    multi_var_est = imfilter(multi_var_est,h,'symmetric');
                case 'average'
                    h = fspecial(EmbeddingParam.Filter.VarPostFilterType,EmbeddingParam.Filter.VarPostFilterSize);
                    multi_var_est = imfilter(multi_var_est,h,'symmetric');
                case 'anisodiff'
                    multi_var_est = anisodiff2D(multi_var_est, 20, 1/7, 0.1, 2);
            end
        end
        FisherInformation = FisherInformation + (1./multi_var_est.^2);
    end
end

if EmbeddingParam.Filter.FIPostFilter
    switch (EmbeddingParam.Filter.FIPostFilterType)
        case 'gaussian'
            h = fspecial(EmbeddingParam.Filter.FIPostFilterType,EmbeddingParam.Filter.FIPostFilterSize,EmbeddingParam.Filter.FIPostFilterSigma);
            FisherInformation = imfilter(FisherInformation,h,'symmetric');
        case 'average'
            h = fspecial(EmbeddingParam.Filter.FIPostFilterType,EmbeddingParam.Filter.FIPostFilterSize);
            FisherInformation = imfilter(FisherInformation,h,'symmetric');
        case 'anisodiff'
            FisherInformation = anisodiff2D(FisherInformation, 20, 1/7, 0.1, 2);
    end
end

I = FisherInformation(:)';

[beta, Cost] = betas_ternary(X,I,alpha);    % Embedding change probabilities
%deflectionCoef = sum(beta.^2 .* I);
Cost = reshape(Cost, size(X));
% profile(X==255) = 10^10;
% profile(X==0) = 10^10;
% Y = STCembedding(X,profile,profile,alpha*numel(X));
if EmbeddingParam.Filter.CostsPostFilter
    switch (EmbeddingParam.Filter.CostsPostFilterType)
        case 'gaussian'
            h = fspecial(EmbeddingParam.Filter.CostsPostFilterType,EmbeddingParam.Filter.CostsPostFilterSize,EmbeddingParam.Filter.CostsPostFilterSigma);
            rho = imfilter(Cost,h,'symmetric');
        case 'average'
            h = fspecial(EmbeddingParam.Filter.CostsPostFilterType,EmbeddingParam.Filter.CostsPostFilterSize);
            rho = imfilter(Cost,h,'symmetric');
        case 'anisodiff'
            rho = anisodiff2D(Cost, 20, 1/7, 0.1, 2);
    end
    
    % adjust embedding costs
    wetCost = 10^10;
    rho(rho > wetCost) = wetCost; % threshold on the costs
    rho(isnan(rho)) = wetCost; % if all xi{} are zero threshold the cost
    rhoP1 = rho;
    rhoM1 = rho;
    rhoP1(X==255) = wetCost; % do not embed +1 if the pixel has max value
    rhoM1(X==0) = wetCost; % do not embed -1 if the pixel has min value
    Cost = rho;
    
    [Y, pChange] = EmbeddingSimulator(X, rhoP1, rhoM1, EmbeddingParam.Payload*numel(X), false);
    change_rate = sum(sum(X~=Y))/numel(X); % Computing the change rate
else
    Y = X;
    beta = 2 * beta;
    r = rand(1,numel(X));
    Modif = (r < beta);                % Cover elements to be modified
    r = rand(1,numel(X));
    Y(Modif) = X(Modif) + 2*(round(r(Modif))) - 1; % Modifying X by +-1
    Y(Y>255) = 254;                    % Taking care of boundary cases
    Y(Y<0)   = 1;
    change_rate = sum(Modif)/numel(X); % Computing the change rate
    pChange = reshape(beta,size(X));
end

end

function Residual = wavden2(X,Nlvl,thrdiv)
[thr,sh,keepapp] = ddencmp('den','wv',X);
Residual = wdencmp('gbl',X,'db8',Nlvl,thr/thrdiv,sh,keepapp);
end

function [y, pChange] = EmbeddingSimulator(x, rhoP1, rhoM1, m, fixEmbeddingChanges)

    n = numel(x);   
    lambda = calc_lambda(rhoP1, rhoM1, m, n);
    pChangeP1 = (exp(-lambda .* rhoP1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));
    pChangeM1 = (exp(-lambda .* rhoM1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));
    pChange = pChangeP1 + pChangeM1; 
    if fixEmbeddingChanges == 1
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',139187));
    else
        RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
    end
    randChange = rand(size(x));
    y = x;
    y(randChange < pChangeP1) = y(randChange < pChangeP1) + 1;
    y(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) = y(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) - 1;
    
    function lambda = calc_lambda(rhoP1, rhoM1, message_length, n)

        l3 = 1e+3;
        m3 = double(message_length + 1);
        iterations = 0;
        while m3 > message_length
            l3 = l3 * 2;
            pP1 = (exp(-l3 .* rhoP1))./(1 + exp(-l3 .* rhoP1) + exp(-l3 .* rhoM1));
            pM1 = (exp(-l3 .* rhoM1))./(1 + exp(-l3 .* rhoP1) + exp(-l3 .* rhoM1));
            m3 = ternary_entropyf(pP1, pM1);
            iterations = iterations + 1;
            if (iterations > 10)
                lambda = l3;
                return;
            end
        end        
        
        l1 = 0; 
        m1 = double(n);        
        lambda = 0;
        
        alpha = double(message_length)/n;
        % limit search to 30 iterations
        % and require that relative payload embedded is roughly within 1/1000 of the required relative payload        
        while  (double(m1-m3)/n > alpha/1000.0 ) && (iterations<30)
            lambda = l1+(l3-l1)/2; 
            pP1 = (exp(-lambda .* rhoP1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));
            pM1 = (exp(-lambda .* rhoM1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));
            m2 = ternary_entropyf(pP1, pM1);
    		if m2 < message_length
    			l3 = lambda;
    			m3 = m2;
            else
    			l1 = lambda;
    			m1 = m2;
            end
    		iterations = iterations + 1;
        end
    end
    
    function Ht = ternary_entropyf(pP1, pM1)
        p0 = 1-pP1-pM1;
        P = [p0(:); pP1(:); pM1(:)];
        H = -((P).*log2(P));
        H((P<eps) | (P > 1-eps)) = 0;
        Ht = sum(H);
    end
end
